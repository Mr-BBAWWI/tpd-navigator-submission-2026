import io
import unittest
import xml.etree.ElementTree as ET

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from packages.science.parent_fingerprint import _svg


def _rdkit_bytes_graph(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise AssertionError(f"failed to parse fixture SMILES: {smiles}")
    for number, atom in enumerate(mol.GetAtoms(), 1):
        atom.SetAtomMapNum(number)
    embedded = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(embedded, randomSeed=19) != 0:
        raise AssertionError("fixture embedding failed")
    source = Chem.RemoveHs(embedded)

    payload = (Chem.MolToMolBlock(source) + "\n$$$$\n").encode("ascii")
    supplier = Chem.ForwardSDMolSupplier(io.BytesIO(payload), removeHs=False)
    parsed = next(iter(supplier), None)
    if parsed is None:
        raise AssertionError("failed to parse RDKit bytes fixture")
    return parsed


def _atom_state(mol: Chem.Mol):
    return tuple(
        (
            atom.GetAtomicNum(),
            atom.GetFormalCharge(),
            atom.GetAtomMapNum(),
            atom.GetIsotope(),
            int(atom.GetChiralTag()),
            atom.GetNoImplicit(),
            atom.GetNumExplicitHs(),
            tuple(
                sorted(
                    (name, atom.GetProp(name))
                    for name in atom.GetPropNames(
                        includePrivate=True,
                        includeComputed=True,
                    )
                )
            ),
        )
        for atom in mol.GetAtoms()
    )


def _svg_text(svg: str) -> tuple[ET.Element, str]:
    root = ET.fromstring(svg)
    return root, "".join(root.itertext())


class ParentChargeSvgTests(unittest.TestCase):
    def test_actual_charged_and_neutral_source_labels_are_distinct(self):
        charged = _rdkit_bytes_graph("C[NH2+]C")
        neutral = _rdkit_bytes_graph("CNC")

        charged_root, charged_text = _svg_text(_svg(charged))
        neutral_root, neutral_text = _svg_text(_svg(neutral))

        self.assertTrue(charged_root.tag.endswith("svg"))
        self.assertTrue(neutral_root.tag.endswith("svg"))
        self.assertIn("N:2 H2 +", charged_text)
        self.assertIn("N:2 H1", neutral_text)
        self.assertNotIn("N:2 H1 +", neutral_text)
        self.assertNotEqual(charged_text, neutral_text)
        self.assertIn(
            "출처 구조의 형식전하·결합 수소 표시, 상태 승인 아님",
            charged_text,
        )
        self.assertIn(
            "출처 구조의 형식전하·결합 수소 표시, 상태 승인 아님",
            neutral_text,
        )

    def test_svg_generation_does_not_change_source_properties_or_coordinates(self):
        source = _rdkit_bytes_graph("C[NH2+]C")
        before_properties = _atom_state(source)
        before_coordinates = np.asarray(
            source.GetConformer().GetPositions(), dtype=float
        ).copy()
        before_is_3d = source.GetConformer().Is3D()

        svg = _svg(source)
        ET.fromstring(svg)

        self.assertEqual(_atom_state(source), before_properties)
        self.assertTrue(
            np.array_equal(
                np.asarray(source.GetConformer().GetPositions(), dtype=float),
                before_coordinates,
            )
        )
        self.assertEqual(source.GetConformer().Is3D(), before_is_3d)
        self.assertTrue(source.GetConformer().Is3D())


if __name__ == "__main__":
    unittest.main()
