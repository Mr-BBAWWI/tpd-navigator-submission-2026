from html.parser import HTMLParser
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from packages.science import parent_fingerprint


def _packet():
    contact = {
        "key": "aromatic_contact_geometry|L:7|A:42:PHE:CG",
        "protein_atom_id": "A:42:PHE:CG",
        "ligand_atom_map": 7,
        "interaction_kind": "aromatic_contact_geometry",
        "heavy_distance_A": 3.61,
        "directional_hydrogen": {
            "distance_HA_A": 2.1,
            "angle_DHA_deg": 151.0,
            "tested_hydrogens": ["H7"],
            "confirmed_geometry": True,
            "ambiguity": "protein donor hydrogen absent",
        },
        "protein_hydrogens_supplied": False,
        "chemical_role_ambiguity": True,
    }
    return {
        "format": "parent-fingerprint/test",
        "structure_svgs": {
            "source_charged_parent": "charged-parent.svg",
            "neutral_design_parent": "neutral-parent.svg",
        },
        "parents": [{
            "parent_role": "source_charged_parent",
            "states": [{
                "state_index": 1,
                "formal_charge": 1,
                "origin": "fixture <unsafe>",
                "contacts": [contact],
                "files": {
                    "raw_sdf": "charged/state01.raw.sdf",
                    "hydrogen_optimized_sdf": "charged/state01.opt.sdf",
                },
            }],
            "common_across_enumerated_states": [contact],
            "proposed_protected_interaction_set": {
                "status": "calculated_proposal_not_expert_approved",
                "contacts": [contact],
            },
        }],
    }


def _sixteen_state_packet():
    packet = _packet()
    contact = packet["parents"][0]["states"][0]["contacts"][0]
    parents = []
    for role, folder in (
        ("source_charged_parent", "charged"),
        ("neutral_design_parent", "neutral"),
    ):
        states = []
        for index in range(1, 9):
            states.append({
                "state_index": index,
                "formal_charge": 1 if role == "source_charged_parent" else 0,
                "origin": "fixture <unsafe>",
                "contacts": [contact],
                "files": {
                    "raw_sdf": f"{folder}/state{index:02d}.raw.sdf",
                    "hydrogen_optimized_sdf": f"{folder}/state{index:02d}.opt.sdf",
                },
            })
        parents.append({
            "parent_role": role,
            "states": states,
            "common_across_enumerated_states": [contact],
            "proposed_protected_interaction_set": {
                "status": "calculated_proposal_not_expert_approved",
                "contacts": [contact],
            },
        })
    packet["parents"] = parents
    return packet


class _EmbeddedScriptParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.current_id = None
        self.parts = {}
        self.script_ids = []

    def handle_starttag(self, tag, attrs):
        if tag == "script":
            self.current_id = dict(attrs).get("id", "__application__")
            self.script_ids.append(self.current_id)
            self.parts.setdefault(self.current_id, [])

    def handle_data(self, data):
        if self.current_id is not None:
            self.parts[self.current_id].append(data)

    def handle_endtag(self, tag):
        if tag == "script":
            self.current_id = None

    def text(self, identifier):
        return "".join(self.parts[identifier])


class ParentPacketHtmlTests(unittest.TestCase):
    def render(self, packet):
        return parent_fingerprint._html_report(packet)

    def parse_scripts(self, html):
        parser = _EmbeddedScriptParser()
        parser.feed(html)
        parser.close()
        return parser

    def test_uses_actual_generated_and_source_paths(self):
        html = self.render(_packet())

        self.assertIn("charged-parent.svg", html)
        self.assertIn("neutral-parent.svg", html)
        self.assertIn("charged/state01.raw.sdf", html)
        self.assertIn("charged/state01.opt.sdf", html)
        self.assertIn("최적화 전 원본 SDF(수소 포함 여부를 단정하지 않음)", html)
        self.assertIn("중원자 고정 및 추가 수소 최적화 SDF", html)
        self.assertIn("생성된 2D 원자 맵 구조", html)

    def test_is_read_only_escaped_and_diagnostic_only(self):
        html = self.render(_packet())

        self.assertNotIn("fixture <unsafe>", html)
        self.assertIn("fixture \\u003cunsafe\\u003e", html)
        self.assertNotIn(".innerHTML", html)
        self.assertNotIn("insertAdjacentHTML", html)
        self.assertNotIn("<button", html)
        self.assertIn("방향족 원자 근접성(π-스태킹 확인이 아님)", html)
        self.assertIn("리간드-H 기하 진단이며 단백질 H는 없음", html)
        self.assertIn("생물학적 승인이 아님", html)
        self.assertIn(
            "전문가 검토 대기 중인 제안일 뿐입니다. 이 패킷은 보호 기준을 확정하지 않습니다.",
            html,
        )
        self.assertIn("common.length", html)

    def test_malicious_json_cannot_end_data_script_or_add_controls(self):
        packet = _packet()
        malicious = "</script><script id='injected'>alert('efficacy')</script><button>approve</button>"
        packet["efficacy_claim"] = malicious
        packet["parents"][0]["states"][0]["origin"] = malicious

        html = self.render(packet)
        parser = self.parse_scripts(html)
        parsed = json.loads(parser.text("data"))

        self.assertEqual(parser.script_ids, ["data", "__application__"])
        self.assertNotIn("injected", parser.parts)
        self.assertEqual(parsed["efficacy_claim"], malicious)
        self.assertEqual(parsed["parents"][0]["states"][0]["origin"], malicious)
        self.assertNotIn("</script><script id='injected'>", html)
        self.assertNotIn("<button>approve</button>", html)
        self.assertIn("\\u003c/script\\u003e", parser.text("data"))

    def test_actual_output_has_no_forbidden_controls_and_parses_payload(self):
        packet = _sixteen_state_packet()
        packet["parents"][0]["states"][0]["origin"] = "한국 경로\x00fixture"
        html = self.render(packet)
        forbidden = [
            (index, ord(character))
            for index, character in enumerate(html)
            if ord(character) < 32 and character not in "\n\t\r"
        ]

        self.assertEqual(forbidden, [])
        parser = self.parse_scripts(html)
        parsed = json.loads(parser.text("data"))
        self.assertEqual(
            sum(len(parent["states"]) for parent in parsed["parents"]),
            16,
        )
        self.assertEqual(
            parsed["structure_svgs"],
            {
                "source_charged_parent": "charged-parent.svg",
                "neutral_design_parent": "neutral-parent.svg",
            },
        )
        self.assertEqual(
            parsed["parents"][0]["states"][0]["origin"],
            "한국 경로\x00fixture",
        )
        self.assertNotIn("\x00", html)
        self.assertGreaterEqual(html.count("charged-parent.svg"), 1)
        self.assertGreaterEqual(html.count("neutral-parent.svg"), 1)
        self.assertIn("charged/state08.opt.sdf", html)
        self.assertIn("neutral/state08.opt.sdf", html)

    def test_rejected_directional_hydrogen_records_render_their_supplied_fields(self):
        node = shutil.which("node")
        if node is None:
            self.skipTest("Node.js is not installed")

        html = self.render(_packet())
        script = self.parse_scripts(html).text("__application__")
        functions = script[
            script.index("function compactValue"):script.index("function drawFigures")
        ]
        row = {
            "protein_atom_id": "A:1421:TYR:OH",
            "directional_hydrogen": {
                "distance_HA_A": None,
                "angle_DHA_deg": None,
                "tested_hydrogens": [{
                    "distance_HA_A": 3.1686,
                    "angle_DHA_deg": 50.212,
                }],
                "confirmed_geometry": False,
                "ambiguity": None,
            },
            "protein_hydrogens_supplied": False,
        }
        program = functions + "\nprocess.stdout.write(JSON.stringify(directionLines(" + json.dumps(row) + ")));"
        result = subprocess.run(
            [node, "-e", program],
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        lines = json.loads(result.stdout)
        self.assertEqual(lines, [
            "검사한 H 1: H···A: 3.1686 Å; D–H···A: 50.212°",
            "방향성 기하가 확인되지 않음",
        ])
        self.assertNotIn("[object Object]", result.stdout)
        self.assertIn("row.protein_atom_id", script)
        self.assertIn("node.textContent=String(text)", script)
        self.assertIn("add(directional,'div',line)", script)

    def test_actual_application_script_is_valid_javascript(self):
        node = shutil.which("node")
        if node is None:
            self.skipTest("Node.js is not installed")

        html = self.render(_sixteen_state_packet())
        parser = self.parse_scripts(html)
        script = parser.text("__application__")

        self.assertNotIn("RegExp(", script)
        self.assertNotIn("/[", script)
        self.assertIn("charCodeAt", script)
        self.assertIn("String.fromCharCode(92)", script)

        with tempfile.TemporaryDirectory(prefix="부모-패킷-") as temporary_directory:
            script_path = Path(temporary_directory) / "보고서-script.js"
            script_path.write_text(script, encoding="utf-8")
            result = subprocess.run(
                [node, "--check", str(script_path)],
                capture_output=True,
                text=True,
                check=False,
            )

        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
