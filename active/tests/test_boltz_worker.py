import json
import tempfile
import unittest
from pathlib import Path

import yaml

from packages.science.boltz_worker import (
    BoltzWorkerError,
    _prediction_sequences,
    assess_outputs,
    build_command,
    discover_outputs,
    prepare_input,
)

import gemmi


SYNTHETIC_CIF = """data_synthetic
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer 'DNA damage-binding protein 1 DDB1 synthetic fixture'
2 polymer 'Cereblon synthetic fixture'
3 polymer 'Bromodomain-containing protein 4 BRD4 synthetic fixture'
4 non-polymer 'RN6 synthetic ligand'
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_seq_one_letter_code_can
1 'polypeptide(L)' AAA
2 'polypeptide(L)' GGG
3 'polypeptide(L)' ACDE
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 2
C 3
E 4
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 1 0 0 0 1 X 1 1
ATOM 2 C CA . ALA B 2 1 0 1 0 1 Y 1 1
ATOM 3 C CA . ALA C 3 1 0 0 1 1 Z 1 1
HETATM 4 C C1 . RN6 E 4 . 1 1 1 1 L 1 1
HETATM 5 O O1 . RN6 E 4 . 2 1 1 1 L 1 1
#
"""


class BoltzWorkerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.reference = self.root / "6BOY.cif"
        self.reference.write_text(SYNTHETIC_CIF, encoding="utf-8")
        self.metadata = self.root / "meta.json"
        self.metadata.write_text(json.dumps({
            "pdb": "6BOY", "ccd": "RN6", "e3_type": "CRBN",
            "reference_target_chain": "C", "reference_e3_chain": "B",
            "reference_ligand_chain": "E",
        }), encoding="utf-8")

    def tearDown(self):
        self.temp.cleanup()

    def test_seed_zero_is_passed_to_confirmed_official_option(self):
        command = build_command(
            executable=Path("boltz.exe"), input_yaml=Path("input.yaml"),
            out_dir=Path("out"), cache=Path("cache"), checkpoint=Path("weights.ckpt"),
            accelerator="gpu", seed=0, msa_mode="single_sequence", max_msa_seqs=32,
            help_has_seed=True,
        )
        index = command.index("--seed")
        self.assertEqual(command[index + 1], "0")
        self.assertNotIn("--use_potentials", command)
        self.assertNotIn("--use_msa_server", command)

    def test_prepared_input_has_no_reference_template_or_coordinate_leakage(self):
        destination = self.root / "input.yaml"
        packet = prepare_input(self.reference, self.metadata, destination, "single_sequence")
        document = yaml.safe_load(destination.read_text(encoding="utf-8"))
        self.assertEqual(set(document), {"version", "sequences"})
        self.assertEqual([next(iter(row)) for row in document["sequences"]],
                         ["protein", "protein", "ligand"])
        self.assertEqual(document["sequences"][2]["ligand"], {"id": "L", "ccd": "RN6"})
        for row in document["sequences"][:2]:
            protein = row["protein"]
            self.assertEqual(protein["msa"], "empty")
            self.assertEqual(set(protein), {"id", "sequence", "msa"})
        serialized = destination.read_text(encoding="utf-8").lower()
        self.assertIn("msa: empty", serialized)
        for forbidden in ("template", "constraint", "coordinate", "6boy.cif", "potential"):
            self.assertNotIn(forbidden, serialized)
        self.assertEqual(packet["explicit_exclusions"]["DDB1"], "입력에서 제외됨")
        self.assertEqual(packet["reference_chain_mapping"]["target"]["auth_asym_id"], "Z")

    def test_server_mode_is_explicit_and_does_not_write_empty_msa(self):
        destination = self.root / "server.yaml"
        prepare_input(self.reference, self.metadata, destination, "server")
        document = yaml.safe_load(destination.read_text(encoding="utf-8"))
        self.assertNotIn("msa", document["sequences"][0]["protein"])
        command = build_command(
            executable=Path("boltz.exe"), input_yaml=destination, out_dir=Path("out"),
            cache=Path("cache"), checkpoint=Path("weights.ckpt"), accelerator="cpu",
            seed=2, msa_mode="server", max_msa_seqs=8, help_has_seed=True,
        )
        self.assertIn("--use_msa_server", command)

    def test_missing_prediction_is_classified_not_success(self):
        output = self.root / "empty_output"
        output.mkdir()
        found = discover_outputs(output)
        self.assertEqual(found["status"], "missing_prediction")
        packet = prepare_input(self.reference, self.metadata, self.root / "missing.yaml",
                               "single_sequence")
        assessed = assess_outputs(output, self.reference, packet)
        self.assertEqual(assessed["status"], "missing_prediction")
        self.assertNotEqual(assessed["status"], "success")

    def test_malformed_confidence_never_becomes_success(self):
        output = self.root / "bad_output"
        output.mkdir()
        (output / "prediction.cif").write_text("not a cif", encoding="utf-8")
        (output / "confidence_prediction.json").write_text("{broken", encoding="utf-8")
        packet = prepare_input(self.reference, self.metadata, self.root / "bad.yaml",
                               "single_sequence")
        assessed = assess_outputs(output, self.reference, packet)
        self.assertEqual(assessed["status"], "inspection_failed")
        self.assertIn("CONFIDENCE_JSON_INVALID", assessed["reason"])

    def test_ambiguous_prediction_files_are_rejected(self):
        output = self.root / "ambiguous"
        output.mkdir()
        (output / "a.cif").write_text("data_a\n", encoding="utf-8")
        (output / "b.cif").write_text("data_b\n", encoding="utf-8")
        (output / "confidence_a.json").write_text('{"score": 0.5}', encoding="utf-8")
        found = discover_outputs(output)
        self.assertEqual(found["status"], "malformed_prediction")
        self.assertEqual(found["reason"], "PREDICTION_CIF_AMBIGUOUS")

    def _prediction_block(self, canonical="XXXX", seq_rows=None):
        if seq_rows is None:
            seq_rows = [(1, "ALA", "n"), (2, "CYS", "n"),
                        (3, "ASP", "n"), (4, "GLU", "n")]
        body = "\n".join(f"1 {number} {monomer} {hetero}"
                         for number, monomer, hetero in seq_rows)
        text = f"""data_prediction
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_seq_one_letter_code_can
1 'polypeptide(L)' {canonical}
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
{body}
loop_
_struct_asym.id
_struct_asym.entity_id
T 1
#
"""
        path = self.root / "prediction-sequence.cif"
        path.write_text(text, encoding="utf-8")
        return gemmi.cif.read_file(str(path)).sole_block()

    def test_prediction_sequence_uses_poly_seq_when_canonical_is_all_x(self):
        self.assertEqual(_prediction_sequences(self._prediction_block()), {"T": "ACDE"})

    def test_prediction_sequence_rejects_conflicting_canonical_value(self):
        with self.assertRaisesRegex(BoltzWorkerError, "PREDICTION_SEQUENCE_SOURCES_DISAGREE"):
            _prediction_sequences(self._prediction_block(canonical="ACDF"))

    def test_prediction_sequence_rejects_missing_poly_seq(self):
        with self.assertRaisesRegex(BoltzWorkerError, "PREDICTION_ENTITY_POLY_SEQ_MISSING"):
            _prediction_sequences(self._prediction_block(seq_rows=[]))

    def test_prediction_sequence_rejects_duplicate_index(self):
        rows = [(1, "ALA", "n"), (1, "CYS", "n")]
        with self.assertRaisesRegex(BoltzWorkerError, "PREDICTION_ENTITY_POLY_SEQ_DUPLICATE_INDEX"):
            _prediction_sequences(self._prediction_block(seq_rows=rows))

    def test_prediction_sequence_rejects_gap_heterogeneous_and_unknown(self):
        cases = [
            ([(1, "ALA", "n"), (3, "ASP", "n")], "PREDICTION_ENTITY_POLY_SEQ_GAP"),
            ([(1, "ALA", "y")], "PREDICTION_ENTITY_POLY_SEQ_HETEROGENEOUS"),
            ([(1, "MSE", "n")], "PREDICTION_ENTITY_POLY_SEQ_UNKNOWN_MONOMER"),
        ]
        for rows_, reason in cases:
            with self.subTest(reason=reason), self.assertRaisesRegex(BoltzWorkerError, reason):
                _prediction_sequences(self._prediction_block(seq_rows=rows_))


if __name__ == "__main__":
    unittest.main()
