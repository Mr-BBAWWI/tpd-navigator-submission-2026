"""Annotated 2D review pictures tied to existing SDF atom maps."""
import io

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from .handoff import check, validate_material
from .molecules import canonical


def attachment_depictions(snapshot, material):
    validate_material(material)
    files = {}
    for compound in material["h1"]["compounds"]:
        name = compound["files"]["sdf_2d"]["path"]
        mols = list(Chem.ForwardSDMolSupplier(io.BytesIO(snapshot.files[name])))
        check(len(mols) == 1 and mols[0] is not None, "DEPICTION_SDF_INVALID")
        mol = mols[0]
        check(canonical(mol) == compound["identity"]["canonical_isomeric_smiles"], "DEPICTION_MOLECULE_MISMATCH")
        atom = compound["attachment_atom"]
        indexes = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum() == atom["atom_map"]]
        check(len(indexes) == 1, "DEPICTION_ATOM_MAP_MISMATCH")
        idx = indexes[0]
        order_index = atom.get("assembled_sdf_index_one_based", atom["rdkit_index_zero_based"] + 1)
        check(idx + 1 == order_index, "DEPICTION_SDF_ORDER_MISMATCH")
        check(mol.HasProp("ccd_atom_names_in_sdf_order")
              and mol.GetProp("ccd_atom_names_in_sdf_order").split(",")[idx] == atom["ccd_atom_id"], "DEPICTION_CCD_ATOM_MISMATCH")
        for a in mol.GetAtoms():
            a.SetAtomMapNum(0)
        mol.GetAtomWithIdx(idx).SetProp("atomNote", atom["ccd_atom_id"])
        rdDepictor.Compute2DCoords(mol)
        drawing = rdMolDraw2D.MolDraw2DSVG(1100, 550)
        drawing.DrawMolecule(mol, highlightAtoms=[idx], highlightAtomColors={idx: (1.0, 0.65, 0.2)},
                             legend=f'{compound["compound_id"]} / {compound["structure"]["ccd"]} / {atom["ccd_atom_id"]} - 2D review only')
        drawing.FinishDrawing()
        files[compound["compound_id"] + "-attachment.svg"] = drawing.GetDrawingText().encode("utf-8")
    return files
