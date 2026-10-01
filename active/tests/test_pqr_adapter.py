import hashlib
import tempfile
import unittest
from pathlib import Path

from packages.science.chemical_states import _actual_pka_entries, read_prepared_protein


def pdb_line(serial, name, residue, chain, sequence, x, y, z, element):
    return (
        f"ATOM  {serial:5d} {name:>4s} {residue:>3s} {chain:1s}{sequence:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
    )


class PqrAdapterTests(unittest.TestCase):
    def test_log_and_captured_stdout_are_actual_pka_sources(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "run.log").write_text(
                "Determinants of pKa\n"
                " ASP1391 A 9.90\n"
                "SUMMARY OF THIS PREDICTION\n"
                " ASP1392 A 4.10\n",
                encoding="utf-8",
            )
            stdout = root / "stdout.bin"
            stdout.write_bytes(
                b"GLU30 B 8.80\n"
                b"SUMMARY OF PKA\n"
                b"GLU31 B 5.20\n"
                b"invalid text\n"
            )
            entries = _actual_pka_entries(root, (stdout,))
            self.assertEqual([
                (item["residue"], item["sequence_number"], item["chain"], item["pKa"])
                for item in entries
            ], [
                ("ASP", 1392, "A", 4.1), ("GLU", 31, "B", 5.2),
            ])
            self.assertNotIn(9.9, [item["pKa"] for item in entries])
            self.assertNotIn(8.8, [item["pKa"] for item in entries])

    def test_prepared_reader_assigns_only_unique_same_residue_hydrogen_parent(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pdb = root / "prepared.pdb"
            pdb.write_text(
                pdb_line(10, "ND1", "HIS", "A", 7, 0.0, 0.0, 0.0, "N")
                + pdb_line(11, "NE2", "HIS", "A", 7, 3.0, 0.0, 0.0, "N")
                + pdb_line(12, "HD1", "HIS", "A", 7, 1.0, 0.0, 0.0, "H"),
                encoding="utf-8",
            )
            digest = hashlib.sha256(pdb.read_bytes()).hexdigest()
            result = read_prepared_protein({
                "output_pdb_path": str(pdb),
                "output_pdb_sha256": digest,
            })
            self.assertEqual(result["protein_hydrogens"][0]["parent_atom_id"], "A:7:HIS:ND1")
            self.assertTrue(result["protein_hydrogens"][0]["parent_topology_reconstructed"])
            nd1 = next(atom for atom in result["protein_atoms"] if atom["label_atom_id"] == "ND1")
            ne2 = next(atom for atom in result["protein_atoms"] if atom["label_atom_id"] == "NE2")
            self.assertEqual(nd1["original_atom_serial"], 10)
            self.assertEqual(nd1["observed_attached_hydrogens"], ["HD1"])
            self.assertFalse(ne2["histidine_other_nitrogen_has_hydrogen"] is False)

    def test_prepared_reader_rejects_changed_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            pdb = Path(temporary) / "prepared.pdb"
            pdb.write_text(pdb_line(1, "CA", "ALA", "A", 1, 0, 0, 0, "C"), encoding="utf-8")
            receipt = {
                "output_pdb_path": str(pdb),
                "output_pdb_sha256": hashlib.sha256(pdb.read_bytes()).hexdigest(),
            }
            pdb.write_text(pdb.read_text(encoding="utf-8") + "REMARK changed\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "provenance hash"):
                read_prepared_protein(receipt)


if __name__ == "__main__":
    unittest.main()
