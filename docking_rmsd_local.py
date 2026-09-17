#!/usr/bin/env python3
"""Symmetry-aware RMSD between a docked multi-model PDBQT and a reference SDF.

GNINA cannot run on this box (its binaries dynamically link libcudart/cusparse/
cublas .so.12) and the free Colab T4 stayed demand-blocked. AutoDock-Vina's pip
wheel (`vina` 1.2.7) is pure CPU and runs locally on 12 cores, but its PDBQT
output is unreadable by RDKit (charge column + AutoDock atom types, e.g. 'A' for
aromatic carbon) and obabel round-trips of it destroy aromaticity.

Key insight: PDBQT columns 1-54 ARE a legal PDB ATOM record - only the trailing
vdW/charge/type fields differ. Truncating to 54 chars yields a block RDKit can
parse with proximityBonding=True, after which AssignBondOrdersFromTemplate
restores correct bond orders and RDKit's GetBestRMS handles molecular symmetry
rigorously.

Two numbers are reported, and they mean different things:
  * receptor-frame RMSD  - NO superposition. This is the honest docking metric:
    does the pose sit where the crystal ligand sits, in the pocket frame?
  * GetBestRMS           - superposition + best-over-symmetry. Measures whether
    the pose *geometry* is right, tolerating a rigid-body offset.

For validation use the receptor-frame number. A pose can have a great GetBestRMS
while sitting in the wrong place; it cannot fake a good receptor-frame RMSD.

Usage: docking_rmsd_local.py <poses.pdbqt> <reference.sdf>
"""
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign

RDLogger.DisableLog("rdApp.*")


def pdbqt_model_to_mol(block, reference):
    """One PDBQT MODEL -> RDKit mol with the reference's bond orders."""
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
        return None  # wrong protonation/element count vs reference
    return Chem.RemoveHs(mol)


def receptor_frame_rmsd(mol, reference):
    """Heavy-atom RMSD in the pocket frame, identity correspondence.

    No superposition, no symmetry search - the strictest possible measure."""
    if mol.GetNumHeavyAtoms() != reference.GetNumHeavyAtoms():
        return None
    pm, rm = mol.GetConformer(), reference.GetConformer()
    p = np.array([[pm.GetAtomPosition(i).x, pm.GetAtomPosition(i).y,
                   pm.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
    r = np.array([[rm.GetAtomPosition(i).x, rm.GetAtomPosition(i).y,
                   rm.GetAtomPosition(i).z] for i in range(reference.GetNumAtoms())])
    d = np.linalg.norm(p - r, axis=1)
    return float(np.sqrt((d ** 2).mean())), float(d.max())


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: docking_rmsd_local.py <poses.pdbqt> <reference.sdf>")
    poses, ref_sdf = Path(sys.argv[1]), Path(sys.argv[2])
    ref = next(m for m in Chem.SDMolSupplier(str(ref_sdf)) if m)
    ref = Chem.RemoveHs(ref)

    text = poses.read_text()
    models = [m for m in text.split("MODEL ")[1:] if m.strip()]
    if not models:
        raise SystemExit(f"no MODEL records in {poses}")
    print(f"reference : {ref_sdf.name} ({ref.GetNumHeavyAtoms()} heavy atoms)")
    print(f"poses     : {poses.name} ({len(models)} models)\n")
    print(f"{'pose':>4} {'receptor-frame':>15} {'GetBestRMS':>11}   verdict")
    best = (np.inf, None)
    for i, block in enumerate(models, 1):
        mol = pdbqt_model_to_mol("MODEL " + block, ref)
        if mol is None:
            print(f"{i:>4}   {'unparsable':>14} {'-':>11}   (element/bond mismatch)")
            continue
        rf = receptor_frame_rmsd(mol, ref)
        try:
            gb = float(rdMolAlign.GetBestRMS(mol, ref))
        except RuntimeError:
            gb = float("nan")
        if rf is None:
            print(f"{i:>4}   {'size mismatch':>14} {gb:>11.3f}")
            continue
        rmsd, mx = rf
        verdict = "PASS (<2.0 A)" if rmsd < 2.0 else "fail"
        print(f"{i:>4}   {rmsd:>10.3f} A {gb:>11.3f}   {verdict}  (max {mx:.2f})")
        if rmsd < best[0]:
            best = (rmsd, i)
    if best[1] is None:
        raise SystemExit("no valid pose could be compared")
    print(f"\nbest receptor-frame RMSD: {best[0]:.3f} A (pose {best[1]})")
    print("criterion <2.0 A = docking reproduces the crystal binding mode")
    return 0 if best[0] < 2.0 else 1


if __name__ == "__main__":
    sys.exit(main())
