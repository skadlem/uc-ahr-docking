#!/usr/bin/env python3
"""Complete Stage A locally: dock the full 424-compound filtered library.

Stage A on Colab covered only chunks 0-3 (300 of 424 compounds); chunks 4-5 were
lost to a corrupt checkpoint and the T4 was then demand-blocked. Local Vina can
now finish the job at the same settings as the original screen (exhaustiveness 2,
no CNN, <=34 heavy atoms) so the results merge cleanly with the Colab data.

This also re-docks everything for a single consistent local dataset, which is
more useful than a partial merge: one scoring function, one receptor prep, one
protonation rule.
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
LIBRARY = ROOT / "screening" / "library_batch1.sdf"
OUTDIR = ROOT / "screening" / "stageA_local"
BOX_CENTER = (160.76, 164.21, 159.80)
BOX_SIZE = 22.0
MAX_HEAVY = 34  # same size filter as the original screen

OBABEL = None


def find_obabel():
    global OBABEL
    if OBABEL:
        return OBABEL
    cand = Path(sys.executable).with_name("obabel")
    if cand.exists():
        OBABEL = str(cand)
        return OBABEL
    import shutil
    OBABEL = shutil.which("obabel")
    return OBABEL


def sdf_to_pdbqt(mol, path):
    tmp = path.with_suffix(".tmp.sdf")
    Chem.MolToMolFile(mol, str(tmp))
    subprocess.run([find_obabel(), str(tmp), "-O", str(path), "-p", "7.4",
                    "--partialcharge", "gasteiger"],
                   capture_output=True, text=True)
    tmp.unlink(missing_ok=True)
    return path.exists() and path.stat().st_size > 0


def strip_salts(smiles):
    """Drop counterions/solvents that obabel writes as disconnected PDBQT fragments.

    A salt like 'C[N+](C)(C)CCO.[Cl-]' becomes two fragments; obabel emits the
    second one as an atom at (0,0,0) after ENDROOT, and Vina rejects the file
    ('Unknown tag in flex residue or ligand'). Keeping only the largest organic
    fragment fixes it. Charges on the parent are preserved."""
    frags = [f for f in smiles.split(".") if f.strip()]
    keep = []
    best = (0, "")
    for f in frags:
        m = Chem.MolFromSmiles(f)
        if m is None:
            continue
        heavy = m.GetNumHeavyAtoms()
        # counterions and common solvents are tiny and/or inorganic
        if heavy <= 3 and not any(a.GetSymbol() == "C" for a in m.GetAtoms()):
            continue
        if heavy > best[0]:
            best = (heavy, f)
        keep.append(f)
    return max(keep, key=len) if keep else smiles


def dock_one(job):
    name, smiles, exh, cpu = job
    pdbqt = OUTDIR / f"{name}.pdbqt"
    try:
        if not pdbqt.exists():
            mol = Chem.MolFromSmiles(strip_salts(smiles))
            if mol is None:
                return name, np.nan, "bad SMILES"
            from rdkit.Chem import AllChem
            mol = Chem.AddHs(mol)
            ps = AllChem.ETKDGv3()
            ps.randomSeed = hash(name) % (2 ** 31)
            AllChem.EmbedMolecule(mol, ps)
            try:
                AllChem.MMFFOptimizeMolecule(mol)
            except Exception:
                pass
            if mol.GetNumHeavyAtoms() > MAX_HEAVY:
                return name, np.nan, "too large"
            if not sdf_to_pdbqt(mol, pdbqt):
                return name, np.nan, "PDBQT failed"
        v = Vina(sf_name="vina", cpu=cpu, verbosity=0)
        v.set_receptor(str(RECEPTOR))
        v.set_ligand_from_file(str(pdbqt))
        v.compute_vina_maps(center=list(BOX_CENTER), box_size=[BOX_SIZE] * 3)
        v.dock(exhaustiveness=exh, n_poses=1)
        e = np.asarray(v.energies())
        score = float((e[:, 0] if e.ndim == 2 else e)[0])
        return name, score, None
    except Exception as exc:
        return name, np.nan, f"{type(exc).__name__}: {exc}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustiveness", type=int, default=2)
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--cpu", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0, help="dock only N (debug)")
    args = ap.parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    jobs = []
    for mol in Chem.SDMolSupplier(str(LIBRARY)):
        if mol is None:
            continue
        cid = mol.GetProp("chembl_id")
        nm = mol.GetProp("name") if mol.HasProp("name") else cid
        if mol.GetNumHeavyAtoms() > MAX_HEAVY:
            continue
        jobs.append((cid, Chem.MolToSmiles(Chem.RemoveHs(mol)),
                     args.exhaustiveness, args.cpu))
    if args.limit:
        jobs = jobs[: args.limit]
    print(f"docking {len(jobs)} compounds (<= {MAX_HEAVY} heavy atoms), "
          f"exh={args.exhaustiveness}, {args.workers} workers x {args.cpu} cpu")
    done = 0
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for name, score, err in pool.map(dock_one, jobs, chunksize=1):
            results.append((name, score, err))
            done += 1
            if done % 50 == 0:
                print(f"  {done}/{len(jobs)}")
    results.sort(key=lambda r: (r[1] if r[1] == r[1] else 1e9))
    csv = OUTDIR / "stageA_local_results.csv"
    with csv.open("w") as fh:
        fh.write("chembl_id,best_dG,error\n")
        for name, score, err in results:
            fh.write(f"{name},{score if score == score else ''},\"{err or ''}\"\n")
    ok = sum(1 for _, s, _ in results if s == s)
    print(f"\nscored {ok}/{len(results)}")
    print(f"{'chembl_id':<14} {'dG':>7}")
    for name, score, err in results[:20]:
        print(f"{name:<14} {score:>7.2f}" + (f"   ERR {err}" if err else ""))
    print(f"\nwrote {csv}")


if __name__ == "__main__":
    main()
