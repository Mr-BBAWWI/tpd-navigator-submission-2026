from __future__ import annotations

import hashlib
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from packages.science.protein_hydrogen_evidence import (
    EvidenceValidationError,
    ProteinHydrogenEvidence,
    prepare_evidence,
)


def pdb_atom(serial, name, residue, chain, sequence, x, y, z, element,
             record="ATOM"):
    return (
        f"{record:<6}{serial:5d} {name:>4s} {residue:>3s} {chain:1s}"
        f"{sequence:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.0:6.2f}{20.0:6.2f}          {element:>2s}\n"
    )


def source_pdb(shift=0.0, omit_ca=False):
    lines = [pdb_atom(1, "N", "ALA", "A", 1, shift, 0, 0, "N")]
    if not omit_ca:
        lines.append(pdb_atom(2, "CA", "ALA", "A", 1, 1.5, 0, 0, "C"))
    lines.append("END\n")
    return "".join(lines)


def source_cif(shift=0.0, omit_ca=False, duplicate_n=False):
    rows = [
        f"ATOM 1 N N ALA A A 1 1 ? {shift:.3f} 0.000 0.000 1.00 . 1"
    ]
    if duplicate_n:
        rows.append(
            f"ATOM 9 N N ALA A A 1 1 ? {shift:.3f} 0.000 0.000 1.00 . 1"
        )
    if not omit_ca:
        rows.append(
            "ATOM 2 C CA ALA A A 1 1 ? 1.500 0.000 0.000 1.00 . 1"
        )
    # Additional deposited protein chain, ligand, and water must remain outside
    # the receptor-PDB equality comparison while staying covered by the CIF hash.
    rows.extend([
        "ATOM 3 N N GLY B B 1 1 ? 8.000 0.000 0.000 1.00 . 1",
        "HETATM 4 C C1 LIG C C 10 10 ? 4.000 3.000 2.000 1.00 . 1",
        "HETATM 5 O O HOH D D 20 20 ? 5.000 3.000 2.000 1.00 . 1",
    ])
    return """data_fixture
_struct.title
;
Deposited fixture with a semicolon-delimited multiline value.
The text is unrelated to atom_site and must parse normally.
;
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.label_alt_id
_atom_site.pdbx_PDB_model_num
""" + "\n".join(rows) + "\n#\n"


def prepared_pdb(shift_n=0.0, omit_ca=False, add_oxt=True,
                 hydrogen_xyz=(-1.0, 0.0, 0.0)):
    lines = [
        pdb_atom(1, "N", "ALA", "A", 1, shift_n, 0, 0, "N"),
        pdb_atom(3, "H", "ALA", "A", 1, *hydrogen_xyz, "H"),
    ]
    if not omit_ca:
        lines.append(pdb_atom(2, "CA", "ALA", "A", 1, 1.5, 0, 0, "C"))
    if add_oxt:
        lines.append(pdb_atom(4, "OXT", "ALA", "A", 1, 2.5, 0, 0, "O"))
    lines.append("END\n")
    return "".join(lines)


def sha(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mocked_run(prepared_text,
               stderr=b"PDB2PQR version 3.6.2\nPROPKA version 3.5.1\n"):
    def implementation(command, cwd, **kwargs):
        output_pdb = next(
            item.split("=", 1)[1]
            for item in command if item.startswith("--pdb-output=")
        )
        output_pqr = Path(command[-1])
        Path(output_pdb).write_text(prepared_text, encoding="utf-8")
        output_pqr.write_text("ATOM synthetic PQR\n", encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout=b"", stderr=stderr)
    return implementation


class ProteinHydrogenEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.pdb = self.root / "input.pdb"
        self.cif = self.root / "input.cif"
        self.executable = self.root / "pdb2pqr"
        self.pdb.write_text(source_pdb(), encoding="utf-8")
        self.cif.write_text(source_cif(), encoding="utf-8")
        self.executable.write_bytes(b"synthetic executable fixture")
        self.executable.chmod(0o755)

    def tearDown(self):
        self.temporary.cleanup()

    def test_semicolon_full_cif_uses_only_matching_protein_subset(self):
        output = self.root / "evidence"
        with mock.patch.object(
            subprocess, "run", side_effect=mocked_run(prepared_pdb())
        ):
            result = prepare_evidence(
                self.pdb, self.cif, output, self.executable, fixed_heavy=True
            )

        self.assertEqual(result["counts"]["source_pdb_heavy_atoms"], 2)
        self.assertEqual(
            result["counts"]["source_cif_selected_matching_protein_heavy_atoms"], 2
        )
        self.assertEqual(
            result["counts"]["source_cif_eligible_heavy_atoms_full_deposition"], 5
        )
        self.assertEqual(result["counts"]["source_heavy_atoms_for_contacts"], 2)
        self.assertEqual(result["counts"]["pending_terminal_oxt_atoms"], 1)
        self.assertTrue(
            result["state_flags"][
                "source_pdb_matches_selected_cif_protein_subset_exactly"
            ]
        )
        self.assertEqual(
            result["source_frame"]["source_cif_scope"]["selected_auth_chains"],
            ["A"],
        )
        self.assertEqual(
            result["source_frame"]["source_cif_scope"][
                "unselected_eligible_heavy_atom_count"
            ],
            3,
        )
        self.assertNotEqual(
            result["sha256"]["source_cif_file"],
            result["sha256"]["source_cif_selected_subset_canonical"],
        )
        self.assertEqual(
            result["sha256"]["source_pdb_selected_subset_canonical"],
            result["sha256"]["source_cif_selected_subset_canonical"],
        )
        self.assertTrue(all(
            atom["label_atom_id"] != "OXT" for atom in result["protein_atoms"]
        ))
        self.assertEqual(
            result["protein_hydrogens"][0]["parent_atom_id"], "A:1:ALA:N"
        )
        self.assertEqual(
            result["protein_hydrogens"][0]["parent_assignment_status"],
            "validated_source_parent",
        )
        self.assertFalse(
            result["state_flags"]["added_hydrogen_orientation_optimized"]
        )
        self.assertTrue(result["state_flags"]["pending_human_review"])
        self.assertFalse(result["state_flags"]["formal_scientific_approval"])
        self.assertEqual(
            result["status"], "computed_diagnostic_pending_human_review"
        )
        self.assertEqual(
            result["fingerprint"]["propka_version_from_actual_output"], "3.5.1"
        )
        self.assertIn("--noopt", result["receipt"]["command"])
        self.assertIn("--nodebump", result["receipt"]["command"])

    def test_shifted_prepared_heavy_atom_is_rejected(self):
        with mock.patch.object(
            subprocess, "run",
            side_effect=mocked_run(prepared_pdb(shift_n=0.001)),
        ):
            with self.assertRaises(EvidenceValidationError) as caught:
                prepare_evidence(
                    self.pdb, self.cif, self.root / "shifted", self.executable
                )
        self.assertIn(
            "source_heavy_atom_position_changed",
            {item["code"] for item in caught.exception.failures},
        )

    def test_missing_prepared_heavy_atom_is_rejected(self):
        with mock.patch.object(
            subprocess, "run",
            side_effect=mocked_run(prepared_pdb(omit_ca=True)),
        ):
            with self.assertRaises(EvidenceValidationError) as caught:
                prepare_evidence(
                    self.pdb, self.cif, self.root / "missing", self.executable
                )
        self.assertIn(
            "source_heavy_atom_missing",
            {item["code"] for item in caught.exception.failures},
        )

    def test_pdb_cif_frame_mismatch_is_rejected(self):
        self.cif.write_text(source_cif(shift=0.001), encoding="utf-8")
        with mock.patch.object(
            subprocess, "run", side_effect=mocked_run(prepared_pdb())
        ):
            with self.assertRaises(EvidenceValidationError) as caught:
                prepare_evidence(
                    self.pdb, self.cif, self.root / "frame", self.executable
                )
        self.assertIn(
            "source_pdb_cif_coordinate_mismatch",
            {item["code"] for item in caught.exception.failures},
        )

    def test_missing_exact_pdb_identity_in_cif_is_rejected(self):
        self.cif.write_text(source_cif(omit_ca=True), encoding="utf-8")
        with mock.patch.object(
            subprocess, "run", side_effect=mocked_run(prepared_pdb())
        ):
            with self.assertRaises(EvidenceValidationError) as caught:
                prepare_evidence(
                    self.pdb, self.cif, self.root / "identity", self.executable
                )
        self.assertIn(
            "source_pdb_atom_missing_from_cif_protein_subset",
            {item["code"] for item in caught.exception.failures},
        )

    def test_duplicate_cif_source_identity_is_rejected(self):
        self.cif.write_text(source_cif(duplicate_n=True), encoding="utf-8")
        with mock.patch.object(
            subprocess, "run", side_effect=mocked_run(prepared_pdb())
        ):
            with self.assertRaises(EvidenceValidationError) as caught:
                prepare_evidence(
                    self.pdb, self.cif, self.root / "duplicate", self.executable
                )
        source_failures = [
            item for item in caught.exception.failures
            if item["code"] == "source_parse_failure"
        ]
        self.assertTrue(source_failures)
        self.assertIn("duplicate_cif_source_identity", source_failures[0]["reason"])

    def test_unknown_hydrogen_parent_is_retained_but_excluded(self):
        text = prepared_pdb(
            add_oxt=False, hydrogen_xyz=(10.0, 10.0, 10.0)
        )
        with mock.patch.object(
            subprocess, "run", side_effect=mocked_run(text)
        ):
            result = prepare_evidence(
                self.pdb, self.cif, self.root / "unknown-parent", self.executable
            )
        hydrogen = result["protein_hydrogens"][0]
        self.assertIsNone(hydrogen["parent_atom_id"])
        self.assertEqual(
            hydrogen["parent_assignment_status"],
            "unknown_excluded_from_directional_tests",
        )
        self.assertTrue(
            result["state_flags"][
                "unknown_hydrogen_parents_excluded_from_directional_tests"
            ]
        )
        kinds = {item.get("kind") for item in result["uncertainties"]}
        self.assertIn("ambiguous_reconstructed_DH_topology", kinds)
        self.assertIn("unassigned_protein_hydrogen", kinds)

    def test_prepared_and_parent_uncertainties_are_both_preserved(self):
        text = prepared_pdb(
            add_oxt=False, hydrogen_xyz=(10.0, 10.0, 10.0)
        )
        with mock.patch.object(
            subprocess, "run", side_effect=mocked_run(text)
        ):
            result = prepare_evidence(
                self.pdb, self.cif, self.root / "uncertainties", self.executable
            )
        kinds = [item.get("kind") for item in result["uncertainties"]]
        self.assertIn("ambiguous_reconstructed_DH_topology", kinds)
        self.assertIn("unassigned_protein_hydrogen", kinds)

    def test_registered_hash_fraud_is_rejected(self):
        output = self.root / "fraud"
        with mock.patch.object(
            subprocess, "run", side_effect=mocked_run(prepared_pdb())
        ):
            prepare_evidence(self.pdb, self.cif, output, self.executable)
        (output / "protein_hydrogen_evidence.json").unlink()
        prepared = output / "pdb2pqr" / "source.prepared.pdb"
        prepared.write_text(
            prepared_pdb() + "REMARK tampered\n", encoding="utf-8"
        )
        input_path = output / "source" / "source.pdb"
        pqr = output / "pdb2pqr" / "source.pqr"
        stdout = output / "pdb2pqr" / "pdb2pqr.stdout.bin"
        stderr = output / "pdb2pqr" / "pdb2pqr.stderr.bin"
        receipt = {
            "input_path": str(input_path),
            "input_sha256": sha(input_path),
            "output_pqr_path": str(pqr),
            "output_sha256": sha(pqr),
            "output_pdb_path": str(prepared),
            "output_pdb_sha256": "0" * 64,
            "stdout_path": str(stdout),
            "stdout_sha256": sha(stdout),
            "stderr_path": str(stderr),
            "stderr_sha256": sha(stderr),
            "executable_path": str(self.executable),
            "executable_sha256": sha(self.executable),
            "preserve_heavy_requested": True,
            "command": [
                str(self.executable), "--with-ph=7.4", "--noopt", "--nodebump"
            ],
            "pKa_entries_from_actual_output": [],
        }
        result = ProteinHydrogenEvidence(
            input_path, output / "source" / "source.cif", receipt, output
        ).evidence()
        self.assertEqual(result["status"], "rejected")
        self.assertIn(
            "receipt_artifact_hash_mismatch",
            {item["code"] for item in result["failures"]},
        )

    def test_optimized_mode_remains_pending_diagnostic(self):
        with mock.patch.object(
            subprocess, "run",
            side_effect=mocked_run(prepared_pdb(add_oxt=False)),
        ):
            result = prepare_evidence(
                self.pdb, self.cif, self.root / "optimized", self.executable,
                fixed_heavy=False,
            )
        self.assertNotIn("--noopt", result["receipt"]["command"])
        self.assertTrue(
            result["state_flags"]["optimized_external_heavy_unchanged"]
        )
        self.assertTrue(result["state_flags"]["computed_diagnostic"])
        self.assertTrue(result["state_flags"]["pending_human_review"])
        self.assertFalse(result["state_flags"]["formal_scientific_approval"])
        self.assertEqual(result["actual_pka_entries"], [])
        self.assertNotIn("population_prediction", result["pKa_population_policy"])


if __name__ == "__main__":
    unittest.main()
