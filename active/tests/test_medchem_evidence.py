from __future__ import annotations

import argparse
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "prepare_medchem_evidence.py"
SPEC = importlib.util.spec_from_file_location("prepare_medchem_evidence", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def fixture_csv(missing: str | None = None) -> bytes:
    """Return a short synthetic CSV fixture; it is not an exact publisher-file copy."""
    header = "Ligands,Smiles,SMARCA2 IC50 (nM),SMARCA4 IC50 (nM),,SMARCA2 DC50 (nM),SMARCA2 Dmax (%),SMARCA4 DC50 (nM),SMARCA4 Dmax (%)"
    smiles = {
        "SMI-6077": "c1ccccc1N1CCCCC1", "SMI-6078": "c1ccccc1N1CCOCC1",
        "SMI-6079": "Clc1ccccc1C1CCCC1", "SMI-6080": "Clc1ccccc1C1CCCCC1",
        "SMI-6084": "Clc1ccccc1C1CCCCCC1", "SMI-6409": "Brc1ccccc1C1CCCC1",
        "SMI-6085": "Brc1ccccc1C1CCCCC1", "SMI-1069": "Clc1ccc(N2CCOCC2)cc1",
        "SMI-1073": "Clc1cccc(C2CCNCC2)c1", "SMI-1074": "Clc1ccc(C2CCNCC2)cc1",
        "SMI-1089": "Brc1ccc(C2CCNCC2)cc1", "SMI-3204": "O=C1Nc2ccccc2C1",
        "SMI-3108": "CCc1ccccc1C(=O)N", "SMI-3100": "Oc1ccc(C2CCCCC2)cc1",
        "SMI-3105": "Oc1nccc(C2CCCCC2)c1", "SMI-3106": "Oc1ccnc(C2CCCCC2)c1",
    }
    rows = []
    for index, compound_id in enumerate(MODULE.PARENT_IDS, 1):
        if compound_id == missing:
            continue
        a = '">10,000"' if compound_id in {"SMI-3105", "SMI-3106"} else str(index * 10)
        b = '">10,000"' if compound_id in {"SMI-3105", "SMI-3106"} else str(index * 10 + 1)
        rows.append(f"{compound_id},{smiles[compound_id]},{a},{b},,,,,")
    rows.append("SMD-6087,CCN1CCC(c2ccc(Cl)cc2)CC1,,, ,5.3,97.2,8.7,97")
    return (header + "\r\r\n" + "\r\r\n".join(rows) + "\r\r\n").encode("cp1252")


def test_actual_cp1252_crcrlf_csv_transformation(tmp_path: Path) -> None:
    source = tmp_path / "jm4c01903_si_006.csv"
    source.write_bytes(fixture_csv())
    records, molecules = MODULE.read_publisher_csv(source)
    assert len(records) == 16
    assert len(molecules) == 16
    assert all(record["record_class"] == "measured_parent_ligand" for record in records)
    record = next(item for item in records if item["compound_id"] == "SMI-3105")
    smarca2, smarca4 = record["measurements"]
    assert smarca2["target"] == "SMARCA2"
    assert smarca4["target"] == "SMARCA4"
    assert smarca4["source_column_label"] == "SMARCA4 IC50 (nM)"
    assert smarca2["raw_string"] == ">10,000"
    assert smarca2["relation"] == ">"
    assert smarca2["value"] == 10000


def test_missing_required_parent_row_fails(tmp_path: Path) -> None:
    source = tmp_path / "source.csv"
    source.write_bytes(fixture_csv(missing="SMI-6080"))
    with unittest.TestCase().assertRaisesRegex(ValueError, "SMI-6080"):
        MODULE.read_publisher_csv(source)


def test_missing_assay_conditions_remain_unknown(tmp_path: Path) -> None:
    source = tmp_path / "source.csv"
    source.write_bytes(fixture_csv())
    records, _ = MODULE.read_publisher_csv(source)
    measurement = records[0]["measurements"][0]
    assert measurement["assay"] == "TR-FRET"
    assert measurement["conditions"] == {"buffer": "unknown", "temperature": "unknown"}
    assert measurement["not_reported_as"] == ["Kd"]


def test_protac_is_not_counted_as_parent(tmp_path: Path) -> None:
    source = tmp_path / "source.csv"
    source.write_bytes(fixture_csv())
    records, _ = MODULE.read_publisher_csv(source)
    assert len(records) == 16
    assert "SMD-6087" not in {record["compound_id"] for record in records}


def test_route_records_are_exact_compound_ids() -> None:
    records = [
        {"compound_id": compound_id, "paper_compound_id": config["paper_compound"]}
        for compound_id, config in MODULE.ROUTE_CONFIG.items()
    ]
    exact_c01 = [record for record in records if record["compound_id"] == "C01"]
    wrong = [record for record in records if record["compound_id"] == "C99"]
    assert len(exact_c01) == 1
    assert wrong == []
    assert exact_c01[0]["paper_compound_id"] == "2"


def test_route_hrms_found_or_obtained_and_same_paragraph_as_nmr() -> None:
    text = """Synthetic methods for the preparation of PROTAC 1 (2)\n\nLong carboxamide name (2)\n\nA mixture of precursor 16 (1 mg) in 0.5M HCl/THF with (1).2HCl, NEt3, NaBH(OAc)3, MgSO4 and DMF was purified by preparative HPLC (5-95 % CH3CN in 0.1 % aq. NH3) to obtain (2) as solid (30.0 mg, 28 % yield).\n\n1H NMR data. HRMS calculated 1.0, found 1.0.\n\nSynthetic methods for the preparation of PROTAC 2 (3)"""
    record = MODULE.route_record("C01", MODULE.ROUTE_CONFIG["C01"], MODULE.paragraphs(text), "fixture", "fixture")
    nmr, hrms = record["characterization"]
    assert nmr["type"] == "NMR"
    assert hrms["type"] == "HRMS"
    assert nmr["paragraph_sha256"] == hrms["paragraph_sha256"]


def test_output_folder_must_be_new(tmp_path: Path) -> None:
    output = tmp_path / "existing"
    output.mkdir()
    args = argparse.Namespace(
        csv=tmp_path / "missing.csv",
        binding_si=tmp_path / "missing-binding.txt",
        route_si=tmp_path / "missing-route.txt",
        route_url="https://example.invalid/source",
        ccd=None,
        output=output,
    )
    with unittest.TestCase().assertRaises(FileExistsError):
        MODULE.build(args)


class MedchemEvidenceTests(unittest.TestCase):
    def _with_tmp_path(self, function) -> None:
        with tempfile.TemporaryDirectory() as directory:
            function(Path(directory))

    def test_actual_cp1252_crcrlf_csv_transformation(self) -> None:
        self._with_tmp_path(test_actual_cp1252_crcrlf_csv_transformation)

    def test_missing_required_parent_row_fails(self) -> None:
        self._with_tmp_path(test_missing_required_parent_row_fails)

    def test_missing_assay_conditions_remain_unknown(self) -> None:
        self._with_tmp_path(test_missing_assay_conditions_remain_unknown)

    def test_protac_is_not_counted_as_parent(self) -> None:
        self._with_tmp_path(test_protac_is_not_counted_as_parent)

    def test_route_records_are_exact_compound_ids(self) -> None:
        test_route_records_are_exact_compound_ids()

    def test_route_hrms_found_or_obtained_and_same_paragraph_as_nmr(self) -> None:
        test_route_hrms_found_or_obtained_and_same_paragraph_as_nmr()

    def test_output_folder_must_be_new(self) -> None:
        self._with_tmp_path(test_output_folder_must_be_new)

    def test_terminal_amine_derivatization_marks_common_nitrogen_on_both_sides(self) -> None:
        left, _ = MODULE.mapped_molecule("CCCN", "left")
        right, _ = MODULE.mapped_molecule("CCCN(C)C", "right")
        records = {
            "left": {"measurements": [], "source": {}},
            "right": {"measurements": [], "source": {}},
        }
        result = MODULE.mcs_pair(
            "left", "right", "terminal_N_derivatization", {"left": left, "right": right}, records
        )
        self.assertGreaterEqual(result["mapping_possibility_count"], 1)
        for possibility in result["mapping_possibilities"]:
            left_n = next(atom.GetAtomMapNum() for atom in left.GetAtoms() if atom.GetSymbol() == "N")
            right_n = next(atom.GetAtomMapNum() for atom in right.GetAtoms() if atom.GetSymbol() == "N")
            self.assertIn(left_n, possibility["left_affected_parent_atom_maps"])
            self.assertIn(right_n, possibility["right_affected_atom_maps"])
            self.assertIn(
                {"left_atom_map": left_n, "right_atom_map": right_n},
                possibility["common_atom_map_correspondences"],
            )


if __name__ == "__main__":
    unittest.main()
