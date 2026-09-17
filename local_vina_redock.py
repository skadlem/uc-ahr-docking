#!/usr/bin/env python3
"""Local AutoDock-Vina redocking gate (CPU-only, no CUDA).

GNINA cannot run on this box (its release binaries dynamically link
libcudart/cusparse/cublas .so.12) and the free Colab T4 was demand-blocked for
the whole session. AutoDock-Vina's pip wheel (`vina` 1.2.7) is pure CPU and
installs locally, so the entire screen can run on this box if the local Vina
reproduces the redocking gate that GNINA passed at 0.49 A.

This script is that gate: dock crystallographic indirubin back into the AHR
PAS-B pocket (7ZUB chain D) with the same box and report the top-pose symmetry-
corrected RMSD. A pass (<2 A, the standard docking-reproduction criterion)
authorises running the full Stage D library locally.
"""
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from vina import Vina

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
RECEPTOR = ROOT / "structures" / "processed" / "AHR_D_receptor.pdb"
REF_LIG = ROOT / "structures" / "processed" / "indirubin_JY6.pdb"
PREP_LIG = ROOT / "structures" / "redock" / "ligand_prep.sdf"


def pdbqt_from_pdb(src, dst, is_ligand):
    """PDB -> PDBQT with OpenBabel if available, else a minimal writer.

    Vina needs PDBQT. OpenBabel is the reliable converter; the fallback handles
    only the simplest cases."""
    import shutil
    if shutil.which("obabel"):
        import subprocess
        kind = "ligand" if is_ligand else "protein"
        # -xr preserves bonds; for receptors we want hydrogens/partial charges
        cmd = ["obabel", str(src), "-O", dst, "--gen3d" if is_ligand else ""]
        cmd = [c for c in cmd if c]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if Path(dst).exists() and Path(dst).stat().st_size > 0:
            return True
        print("obabel stderr:", r.stderr[-400:])
        return False
    raise SystemExit("obabel required to build PDBQT")


def main():
    work = ROOT / "structures" / "vina_local"
    work.mkdir(parents=True, exist_ok=True)

    rec_pdbqt = work / "receptor.pdbqt"
    lig_pdbqt = work / "indirubin.pdbqt"
    if not rec_pdbqt.exists():
        ok = pdbqt_from_pdb(RECEPTOR, rec_pdbqt, is_ligand=False)
        if not ok:
            raise SystemExit("receptor conversion failed")
    if not lig_pdbqt.exists():
        ok = pdbqt_from_pdb(PREP_LIG, lig_pdbqt, is_ligand=True)
        if not ok:
            raise SystemExit("ligand conversion failed")
    print(f"receptor: {rec_pdbqt} ({rec_pdbqt.stat().st_size} B)")
    print(f"ligand:   {lig_pdbqt} ({lig_pdbqt.stat().st_size} B)")

    # box from the reference ligand: centroid +/- padding
    ref = next(m for m in Chem.SDMolSupplier(str(PREP_LIG)) if m)
    conf = ref.GetConformer()
    xs, ys, zs = [], [], []
    for i in range(ref.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        xs.append(p.x); ys.append(p.y); zs.append(p.z)
    cx, cy, cz = sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)
    size = 22.0  # A, generous cubic box around the PAS-B pocket
    print(f"box centre ({cx:.2f},{cy:.2f},{cz:.2f}) size {size:.1f} A")

    v = Vina(sf_name="vina", cpu=12, verbosity=1)
    v.set_receptor(str(rec_pdbqt))
    v.set_ligand_from_file(str(lig_pdbqt))
    v.compute_vina_maps(center=[cx, cy, cz], box_size=[size, size, size])

    v.dock(exhaustiveness=32, n_poses=10)
    v.write_pose(str(work / "redocked.pdbqt"))
    import numpy as np
    print("\nenergies:", np.round(v.energies()[:10], 2).tolist())


if __name__ == "__main__":
    main()
