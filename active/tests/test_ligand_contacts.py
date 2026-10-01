"""Synthetic audit/diagnostic edge cases, no molecular execution or approval."""
import copy
from pathlib import Path
import tempfile
import unittest

from packages.science.evidence_common import file_index, seal, write_new_directory
from packages.science.handoff import encoded, sha
from packages.science.ligand_contacts import analyze_atom_audit, archive_files, project_contacts


def fixture():
    ids = [["L", 1, "", "SYN", n, ""] for n in ("N1", "C2", "H1")] + [["P", 2, "", "ALA", "O", ""]]
    elements = ("N", "C", "H", "O")
    coords = ([0., 0., 0.], [1.4, 0., 0.], [-.9, 0., 0.], [-2.0, 0., 0.])
    rows = [{"atom": key, "element": e, "xyz_A": list(x)} for key, e, x in zip(ids, elements, coords)]
    audit_rows = []
    for i, row in enumerate(rows):
        audit_rows.append({**copy.deepcopy(row), "neighbors": [[ids[j] for j in (1, 2)], [ids[0]], [ids[0]], []][i],
                           "energy_type": ["NH", "CH3", "HNH", "O"][i], "restraint_h_bond_type": ["D", "N", "H", "A"][i],
                           "typing_error": None, "probe_radius_A": 1.05 if i == 2 else 1.55, "probe_acceptor": i == 3,
                           "probe_donor": i == 2, "probe_charge": 0})
    context = {"digest": "1" * 64, "prior_quality_digest": "2" * 64, "data_mode": "synthetic_test",
               "binding": {"compound_id": "synthetic", "molecule_id": "synthetic", "run_id": "synthetic", "model_rank": 0, "input_sha256": "3"*64, "structure_sha256": "4"*64},
               "ligand_prefix": ids[0][:4], "dictionary_atoms": {n: {} for n in ("N1", "C2", "H1")},
               "dictionary_bonds": [{"atoms": ["N1", "C2"]}, {"atoms": ["N1", "H1"]}],
               "heavy_atom_roles": [{"name": "N1", "atom_map": 1, "element": "N", "formal_charge": 0, "rdkit_acceptor": False, "rdkit_donor": True}],
               "protein_name_aliases": [],
               "method": "RDKit H directions, exact dictionary X-ray/neutron H lengths, unchanged protein H and all heavy coordinates",
               "models": {m: {"sha256": sha(m.encode()), "atoms": copy.deepcopy(rows), "target_chains": ["P"], "h_bond_lengths": [{}]} for m in ("xray", "nuclear")}}
    files = {}
    for mode in ("xray", "nuclear"):
        for direction in ("forward", "reverse"):
            prefix = f"runs/{mode}/{direction}"
            source, target = (["L"], ["P"]) if direction == "forward" else (["P"], ["L"])
            audit = {"format": "tpd-probe-atom-audit/0.1.0-draft", "cctbx_version": "2025.11", "input_sha256": context["models"][mode]["sha256"],
                     "hydrogen_convention": mode, "source_chains": source, "target_chains": target, "atoms": copy.deepcopy(audit_rows)}
            s, t = (ids[2], ids[3]) if direction == "forward" else (ids[3], ids[2])
            atom = lambda v: dict(zip(("chainID", "resID", "iCode", "resName", "atomName", "alt"), v))
            output = {"flat_results": [{"src": atom(s), "target": atom(t), "group": "1->2", "type": "bo", "gap": -0.5, "dotCount": 2}]}
            files[prefix + "/probe.json"] = encoded(output)
            files[prefix + "/atom-audit.json"] = encoded(audit)
            files[prefix + "/probe.log"] = b"synthetic fixture; no tool execution"
            record = {"exit_code": 0, "error": None, "elapsed_seconds": 0.1, "mode": mode, "source_chains": source, "target_chains": target,
                      "input_sha256": context["models"][mode]["sha256"], "output_sha256": sha(files[prefix + "/probe.json"]),
                      "audit_sha256": sha(files[prefix + "/atom-audit.json"]), "log_sha256": sha(files[prefix + "/probe.log"])}
            files[prefix + "/execution.json"] = encoded(record)
    return context, files, {**audit, "atoms": copy.deepcopy(audit_rows)}


class ContactAuditTests(unittest.TestCase):
    def setUp(self):
        self.context, self.files, self.audit = fixture()
        self.audit.update(hydrogen_convention="xray", input_sha256=self.context["models"]["xray"]["sha256"])

    def test_bidirectional_pairs_are_deduplicated_without_approval(self):
        report = project_contacts(self.context, self.files)
        self.assertEqual(report["status"], "completed_with_limits")
        self.assertEqual([m["bad_overlap_pair_count"] for m in report["models"]], [1, 1])
        self.assertFalse(report["automatic_acceptance"])
        self.assertFalse(report["whole_structure_validated"])
        self.assertTrue(all(m["coverage"] == "incomplete" for m in report["models"]))

    def test_actual_h_donor_is_compared_at_its_parent(self):
        result = analyze_atom_audit(self.context, "xray", self.audit)
        self.assertEqual(result["role_disagreements"], [])
        self.assertEqual(result["status"], "review_required")

    def test_acceptor_disagreement_is_exposed_not_corrected(self):
        self.audit["atoms"][0]["probe_acceptor"] = True
        result = analyze_atom_audit(self.context, "xray", self.audit)
        self.assertEqual(result["role_disagreements"][0]["name"], "N1")
        self.assertTrue(self.audit["atoms"][0]["probe_acceptor"])

    def test_unknown_type_is_incomplete(self):
        self.audit["atoms"][0]["restraint_h_bond_type"] = "?"
        result = analyze_atom_audit(self.context, "xray", self.audit)
        self.assertEqual(result["status"], "incomplete")
        self.assertEqual(len(result["missing_atom_types"]), 1)

    def test_missing_h_bond_rejects_even_successful_tool(self):
        self.audit["atoms"][0]["neighbors"].pop()
        with self.assertRaisesRegex(ValueError, "CONTACT_PROBE_BOND_GRAPH"):
            analyze_atom_audit(self.context, "xray", self.audit)

    def test_unexpected_covalent_link_rejected(self):
        self.audit["atoms"][0]["neighbors"].append(self.audit["atoms"][3]["atom"])
        with self.assertRaisesRegex(ValueError, "CONTACT_UNEXPECTED_COVALENT_LINK"):
            analyze_atom_audit(self.context, "xray", self.audit)

    def test_internal_three_decimal_rounding_is_measured(self):
        self.context["models"]["xray"]["atoms"][0]["xyz_A"][0] = 0.00049
        result = analyze_atom_audit(self.context, "xray", self.audit)
        self.assertAlmostEqual(result["probe_coordinate_rounding_max_A"], 0.00049)

    def test_coordinate_change_beyond_rounding_rejected(self):
        self.audit["atoms"][0]["xyz_A"][0] = 0.001
        with self.assertRaisesRegex(ValueError, "CONTACT_PROBE_CHANGED_MODEL"):
            analyze_atom_audit(self.context, "xray", self.audit)

    def test_coordinate_change_not_at_serialization_precision_rejected(self):
        self.audit["atoms"][0]["xyz_A"][0] = 0.0002
        with self.assertRaisesRegex(ValueError, "CONTACT_PROBE_CHANGED_MODEL"):
            analyze_atom_audit(self.context, "xray", self.audit)

    def test_different_model_binding_rejected(self):
        self.audit["input_sha256"] = "f" * 64
        with self.assertRaisesRegex(ValueError, "CONTACT_AUDIT_BINDING"):
            analyze_atom_audit(self.context, "xray", self.audit)

    def test_missing_atom_rejected(self):
        self.audit["atoms"].pop()
        with self.assertRaisesRegex(ValueError, "CONTACT_AUDIT_ATOM_SET"):
            analyze_atom_audit(self.context, "xray", self.audit)

    def test_charge_mismatch_rejected(self):
        self.audit["atoms"][0]["probe_charge"] = 1
        with self.assertRaisesRegex(ValueError, "CONTACT_PROBE_CHARGE_MISMATCH"):
            analyze_atom_audit(self.context, "xray", self.audit)

    def test_failed_tools_have_no_successful_zero_clash_count(self):
        import json
        for name in list(self.files):
            if name.endswith("execution.json"):
                record = json.loads(self.files[name]); record["exit_code"] = 124; record["error"] = "CONTACT_TOOL_TIMEOUT"
                self.files[name] = encoded(record)
        report = project_contacts(self.context, self.files)
        self.assertEqual(report["status"], "failed_or_partial")
        self.assertTrue(all(m["bad_overlap_pair_count"] is None for m in report["models"]))

    def test_changed_probe_output_is_not_consumed(self):
        self.files["runs/xray/forward/probe.json"] = encoded({"flat_results": []})
        report = project_contacts(self.context, self.files)
        self.assertEqual(report["status"], "failed_or_partial")
        self.assertIn("OUTPUT_HASH", report["models"][0]["directions"]["forward"]["error"])

    def test_archive_rejects_tampering(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "archive"; files = {"data.json": b"{}"}
            write_new_directory(root, {**files, "manifest.json": encoded(seal({"format": "synthetic", "files": file_index(files)}))})
            archive_files(root, "synthetic")
            (root / "data.json").write_bytes(b"modified")
            with self.assertRaisesRegex(ValueError, "CONTACT_ARCHIVE_HASH"):
                archive_files(root, "synthetic")
