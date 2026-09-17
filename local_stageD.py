#!/usr/bin/env python3
"""Local Stage D: dock the 20-compound validation library into AHR PAS-B.

Replaces the Colab T4 run (scripts/cell_stageD.py) now that AutoDock-Vina runs
locally on CPU. Same three-group design:
  BENCH - published AHR modulators (Mosa et al. 2024) for external calibration
  ORTH  - candidates found by the local shape/pharmacophore screens
  CTRL  - indirubin, the 0.51 A redocking positive control

Ligand prep mirrors the receptor: SDF -> obabel -> rigid PDBQT at pH 7.4 with
Gasteiger charges. Obabel's --gen3d is NOT used (the library already has 3D);
obabel is only re-writing the existing geometry into PDBQT form.

Usage: local_stageD.py [--exhaustiveness 16] [--workers 6]
"""
import argparse
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from vina import Vina

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
RECEPTOR = ROOT / "structures" / "vina_local" / "receptor.pdbqt"
LIBRARY = ROOT / "screening" / "validation_library.sdf"
OUTDIR = ROOT / "screening" / "stageD_local"

BOX_CENTER = (160.76, 164.21, 159.80)
BOX_SIZE = 22.0

# subset of Vina output modes that stay near the pocket
N_POSES = 3


OBABEL = None


def find_obabel():
    """Locate obabel next to the active python, else on PATH."""
    global OBABEL
    if OBABEL:
        return OBABEL
    cand = Path(sys.executable).with_name("obabel")
    if cand.exists():
        OBABEL = str(cand)
        return OBABEL
    import shutil
    found = shutil.which("obabel")
    if found:
        OBABEL = found
        return OBABEL
    raise SystemExit("obabel not found; install openbabel-wheel into the venv")


def sdf_to_pdbqt(mol, path):
    """RDKit mol (3D) -> flexible PDBQT with Gasteiger charges.

    Ligands must be written as FLEXIBLE PDBQT (with ROOT/BRANCH torsion records)
    - a rigid -xr conversion emits AutoDock 'A' aromatic-carbon types, which
    Vina's ligand parser rejects ('Unknown tag in flex residue or ligand').
    Geometry is preserved; obabel only adds pH 7.4 hydrogens and charges."""
    tmp = path.with_suffix(".tmp.sdf")
    Chem.MolToMolFile(mol, str(tmp))
    subprocess.run(
        [find_obabel(), str(tmp), "-O", str(path), "-p", "7.4",
         "--partialcharge", "gasteiger"],
        capture_output=True, text=True)
    tmp.unlink(missing_ok=True)
    return path.exists() and path.stat().st_size > 0


def dock_one(args):
    name, group, pdbqt_path, exh, cpu = args
    try:
        v = Vina(sf_name="vina", cpu=cpu, verbosity=0)
        v.set_receptor(str(RECEPTOR))
        v.set_ligand_from_file(str(pdbqt_path))
        v.compute_vina_maps(center=list(BOX_CENTER),
                            box_size=[BOX_SIZE] * 3)
        v.dock(exhaustiveness=exh, n_poses=N_POSES)
        out = OUTDIR / f"{name}.pdbqt"
        if out.exists():
            out.unlink()
        # write_pose() omits the per-pose VINA RESULT remarks for flexible
        # ligands in this wheel; poses() returns the full multi-model text.
        out.write_text(v.poses(n_poses=N_POSES))
        # energies() is (n_poses x n_terms): col 0 = inter+intra, i.e. the
        # reported docking score; other columns are its components.
        e = np.asarray(v.energies())
        if e.ndim == 2:
            scores = [float(x) for x in e[:, 0]]
        else:
            scores = [float(x) for x in e]
        return name, group, scores[0], [round(s, 2) for s in scores[:N_POSES]], None
    except Exception as exc:  # keep the pool alive on a single failure
        return name, group, float("nan"), [], f"{type(exc).__name__}: {exc}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustiveness", type=int, default=16)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--cpu", type=int, default=2)
    args = ap.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    if not RECEPTOR.exists():
        raise SystemExit(f"missing receptor {RECEPTOR}; run scripts/local_vina_redock.py")

    suppl = Chem.SDMolSupplier(str(LIBRARY))
    jobs = []
    for mol in suppl:
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else "x"
        group = mol.GetProp("group") if mol.HasProp("group") else "?"
        pdbqt = OUTDIR / f"{name}.ligand.pdbqt"
        if not pdbqt.exists():
            if not sdf_to_pdbqt(mol, pdbqt):
                print(f"SKIP {name}: PDBQT conversion failed")
                continue
        jobs.append((name, group, pdbqt, args.exhaustiveness, args.cpu))

    print(f"docking {len(jobs)} ligands: exh={args.exhaustiveness}, "
          f"{args.workers} workers x {args.cpu} cpus, {N_POSES} poses each\n")
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(dock_one, jobs, chunksize=1))

    results.sort(key=lambda r: (r[2] if r[2] == r[2] else 1e9))
    csv = OUTDIR / "stageD_local_results.csv"
    with csv.open("w") as fh:
        fh.write("name,group,best_kcal_per_mol,all_poses,error\n")
        for name, group, best, poses, err in results:
            fh.write(f"{name},{group},{best:.2f},\"{poses}\",\"{err or ''}\"\n")
    print(f"\n{'drug':<22} {'group':<7} {'dG (kcal/mol)':>14}")
    for name, group, best, poses, err in results:
        flag = "  <-- ERROR" if err else ""
        print(f"{name:<22} {group:<7} {best:>14.2f}{flag}")
    print(f"\nwrote {csv}")


if __name__ == "__main__":
    main()
