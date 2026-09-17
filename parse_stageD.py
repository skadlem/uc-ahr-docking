#!/usr/bin/env python3
"""Parse Vina pose PDBQT files into a results CSV.

Vina writes `REMARK VINA RESULT:  <affinity> <rmsd_lb> <rmsd_ub>` per pose, so
scores can be recovered from the pose files without re-docking. Also reports the
receptor-frame RMSD of the best pose against the indirubin reference for CTRL.
"""
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / "screening" / "stageD_local"
REF = ROOT / "structures" / "redock" / "ligand_prep.sdf"
NAMES = ROOT / "screening" / "validation_library.sdf"


def chembl_names():
    """CHEMBL id -> (drug name, group) from the library SDF."""
    out = {}
    suppl = Chem.SDMolSupplier(str(NAMES))
    for mol in suppl:
        if mol is None:
            continue
        cid = mol.GetProp("_Name")
        group = mol.GetProp("group") if mol.HasProp("group") else "?"
        out[cid] = (mol.GetProp("drug_name") if mol.HasProp("drug_name") else cid,
                    group)
    return out


def parse_poses(path):
    """[(affinity, rmsd_lb, rmsd_ub)] from a Vina pose PDBQT."""
    scores = []
    for line in path.read_text().split("\n"):
        if line.startswith("REMARK VINA RESULT:"):
            parts = line.split()[3:6]
            scores.append(tuple(float(x) for x in parts))
    return scores


def receptor_frame_rmsd(pose_mol, reference):
    if pose_mol.GetNumHeavyAtoms() != reference.GetNumHeavyAtoms():
        return None
    pm, rm = pose_mol.GetConformer(), reference.GetConformer()
    p = np.array([[pm.GetAtomPosition(i).x, pm.GetAtomPosition(i).y,
                   pm.GetAtomPosition(i).z] for i in range(pose_mol.GetNumAtoms())])
    r = np.array([[rm.GetAtomPosition(i).x, rm.GetAtomPosition(i).y,
                   rm.GetAtomPosition(i).z] for i in range(reference.GetNumAtoms())])
    return float(np.sqrt((np.linalg.norm(p - r, axis=1) ** 2).mean()))


def pose_to_mol(block, reference):
    lines = [l[:54] for l in block.split("\n") if l.startswith(("ATOM", "HETATM"))]
    if not lines:
        return None
    mol = Chem.MolFromPDBBlock("\n".join(lines) + "\nEND\n",
                              sanitize=False, removeHs=True, proximityBonding=True)
    if mol is None:
        return None
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(reference, mol)
    except (RuntimeError, ValueError):
        return None
    return Chem.RemoveHs(mol)


def main():
    ref = next(m for m in Chem.SDMolSupplier(str(REF)) if m)
    ref_h = Chem.RemoveHs(ref)
    names = chembl_names()

    rows = []
    for path in sorted(OUTDIR.glob("*.pdbqt")):
        if path.name.endswith(".ligand.pdbqt"):
            continue
        stem = path.stem
        scores = parse_poses(path)
        if not scores:
            continue
        drug, group = names.get(stem, (stem, "?"))
        # receptor-frame RMSD of the best pose, where reference applies
        rf = None
        try:
            txt = path.read_text()
            block = "MODEL " + txt.split("MODEL ")[1] if "MODEL " in txt else txt
            mol = pose_to_mol(block, ref_h)
            if mol is not None:
                rf = receptor_frame_rmsd(mol, ref_h)
        except Exception:
            pass
        rows.append((drug, group, scores[0][0], scores[0][1], rf,
                     [round(s[0], 2) for s in scores]))

    rows.sort(key=lambda r: r[2])
    out = OUTDIR / "stageD_results.csv"
    with out.open("w") as fh:
        fh.write("drug,group,best_dG,rmsd_lb_vs_best,receptor_frame_rmsd_vs_indirubin,all_pose_dG\n")
        for drug, group, dg, lb, rf, allp in rows:
            rf_s = f"{rf:.3f}" if rf is not None else "n/a"
            fh.write(f"{drug},{group},{dg:.2f},{lb:.2f},{rf_s},\"{allp}\"\n")

    print(f"{'drug':<20} {'group':<6} {'dG':>7} {'RMSD_lb':>8} {'vs indirubin':>13}")
    for drug, group, dg, lb, rf, _ in rows:
        rf_s = f"{rf:.3f}" if rf is not None else "n/a"
        print(f"{drug:<20} {group:<6} {dg:>7.2f} {lb:>8.2f} {rf_s:>13}")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
