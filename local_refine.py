#!/usr/bin/env python3
"""Local refinement of the consensus shortlist at high exhaustiveness.

Stage D established that local Vina reproduces published AHR SAR (indirubin
control 0.428 A). This runs the consensus shortlist at exhaustiveness=64 - the
final ranked poses - and also re-docks the Stage A top hits for a consistent
comparison, since Stage A used GNINA/Vina on Colab.

Also writes per-residue pocket contacts for the best pose of each ligand, which
is the Stage C analysis that Colab never delivered.
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
REF_LIG = ROOT / "structures" / "redock" / "ligand_prep.sdf"
POCKET_PDB = ROOT / "structures" / "processed" / "AHR_D_receptor.pdb"
OUTDIR = ROOT / "screening" / "refine_local"
BOX_CENTER = (160.76, 164.21, 159.80)
BOX_SIZE = 22.0
CONTACT_CUTOFF = 4.5  # heavy-atom distance for a pocket contact

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
    if OBABEL:
        return OBABEL
    raise SystemExit("obabel not found")


def sdf_to_pdbqt(mol, path):
    tmp = path.with_suffix(".tmp.sdf")
    Chem.MolToMolFile(mol, str(tmp))
    subprocess.run([find_obabel(), str(tmp), "-O", str(path), "-p", "7.4",
                    "--partialcharge", "gasteiger"],
                   capture_output=True, text=True)
    tmp.unlink(missing_ok=True)
    return path.exists() and path.stat().st_size > 0


def pocket_residues():
    """Residue name -> nearest-atom position, from the processed receptor PDB."""
    residues = {}
    for line in POCKET_PDB.read_text().split("\n"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        res = line[17:20].strip() + line[21:26].strip()
        if res in ("HOH", "WAT"):
            continue
        pos = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        residues.setdefault(res, []).append(pos)
    return {k: np.array(v) for k, v in residues.items()}


def pose_block(path):
    text = path.read_text()
    if "MODEL " in text:
        return "MODEL " + text.split("MODEL ")[1]
    return text


def pose_contacts(pdbqt_path, residues):
    """Heavy atoms of the best pose -> contacting residues within the cutoff."""
    block = pose_block(pdbqt_path)
    lig_atoms = []
    for line in block.split("\n"):
        if line.startswith(("ATOM", "HETATM")):
            atype = line[77:79].strip().upper()
            if atype in ("HD", "HS", "H"):
                continue
            lig_atoms.append([float(line[30:38]), float(line[38:46]),
                              float(line[46:54])])
    if not lig_atoms:
        return {}
    lig = np.array(lig_atoms)
    contacts = {}
    for res, atoms in residues.items():
        d = np.linalg.norm(lig[:, None, :] - atoms[None, :, :], axis=2).min()
        if d <= CONTACT_CUTOFF:
            contacts[res] = round(float(d), 2)
    return dict(sorted(contacts.items(), key=lambda kv: kv[1]))


def dock_one(job):
    name, sdf_path, exh, cpu = job
    try:
        pdbqt = OUTDIR / f"{name}.pdbqt"
        # name is '<readable>__<chembl_id>'; the ligand must be selected from the
        # (multi-molecule) SDF by that id, not taken as the first entry.
        cid = name.split("__")[1]
        if not pdbqt.exists():
            mol = None
            for m in Chem.SDMolSupplier(str(sdf_path)):
                if m is None:
                    continue
                got = m.GetProp("_Name") if m.HasProp("_Name") else ""
                if got == cid:
                    mol = m
                    break
                got = m.GetProp("name") if m.HasProp("name") else ""
                if got.upper() == cid.upper():
                    mol = m
                    break
            if mol is None:
                return name, np.nan, None, f"{cid} not found in {sdf_path.name}"
            if not sdf_to_pdbqt(mol, pdbqt):
                return name, np.nan, None, "PDBQT conversion failed"
        v = Vina(sf_name="vina", cpu=cpu, verbosity=0)
        v.set_receptor(str(RECEPTOR))
        v.set_ligand_from_file(str(pdbqt))
        v.compute_vina_maps(center=list(BOX_CENTER), box_size=[BOX_SIZE] * 3)
        v.dock(exhaustiveness=exh, n_poses=3)
        out = OUTDIR / f"{name}.pdbqt"
        if out.exists():
            out.unlink()
        out.write_text(v.poses(n_poses=3))
        e = np.asarray(v.energies())
        scores = [float(x) for x in (e[:, 0] if e.ndim == 2 else e)]
        return name, scores[0], scores, None
    except Exception as exc:
        return name, np.nan, None, f"{type(exc).__name__}: {exc}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustiveness", type=int, default=64)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--cpu", type=int, default=2)
    args = ap.parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # consensus shortlist + Stage A leaders, all taken from the libraries on disk
    wanted = {
        "berberine": "screening/library_batch1.sdf",
        "trioxsalen": "screening/library_batch1.sdf",
        "opicapone": "screening/library_batch1.sdf",
        "ramosetron": "screening/validation_library.sdf",
        "triamterene": "screening/validation_library.sdf",
        "methylene_blue": "screening/validation_library.sdf",
        "indalpine": "screening/validation_library.sdf",
    }
    # map the wanted names to the actual ChEMBL ids present in each file.
    # library_batch1 uses the property 'name' (uppercase, e.g. 'BERBERINE');
    # validation_library uses '_Name' with a lowercase chembl id. Match both.
    id_by_name = {}
    for path in set(wanted.values()):
        for mol in Chem.SDMolSupplier(str(ROOT / path)):
            if mol is None:
                continue
            nm = None
            for prop in ("name", "drug_name", "_Name"):
                if mol.HasProp(prop):
                    nm = mol.GetProp(prop).lower().replace(" ", "_")
                    break
            if nm and nm in wanted:
                id_by_name.setdefault(nm, (mol.GetProp("_Name"), path))
    id_by_name["indirubin_ref"] = ("INDIRUBIN_REF", "screening/validation_library.sdf")
    jobs = []
    seen = set()
    for nm, (cid, path) in id_by_name.items():
        key = f"{nm}__{cid}"
        if key in seen:
            continue
        seen.add(key)
        jobs.append((key, ROOT / path, args.exhaustiveness, args.cpu))
    missing = [w for w in wanted if w not in id_by_name]
    if missing:
        print(f"WARNING: not found in libraries: {missing}")
    print(f"refining {len(jobs)} compounds at exh={args.exhaustiveness}, "
          f"{args.workers} workers x {args.cpu} cpus\n")

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(dock_one, jobs, chunksize=1))
    results.sort(key=lambda r: (r[1] if r[1] == r[1] else 1e9))

    residues = pocket_residues()
    csv = OUTDIR / "refine_results.csv"
    with csv.open("w") as fh:
        fh.write("compound,best_dG,all_pose_dG,contacts,error\n")
        for name, best, scores, err in results:
            tag = name.split("__")[0]
            pose = OUTDIR / f"{name}.pdbqt"
            contacts = pose_contacts(pose, residues) if pose.exists() else {}
            fh.write(f"{tag},{best:.2f},\"{[round(s,2) for s in scores] if scores else []}\","
                     f"\"{';'.join(f'{k}:{v}' for k,v in list(contacts.items())[:12])}\","
                     f"\"{err or ''}\"\n")
    print(f"{'compound':<16} {'dG':>7}  nearest pocket contacts")
    for name, best, scores, err in results:
        tag = name.split("__")[0]
        pose = OUTDIR / f"{name}.pdbqt"
        contacts = pose_contacts(pose, residues) if pose.exists() else {}
        cs = ", ".join(f"{k}:{v}" for k, v in list(contacts.items())[:6])
        flag = "  ERROR: " + err if err else ""
        print(f"{tag:<16} {best:>7.2f}  {cs}{flag}")
    print(f"\nwrote {csv}")


if __name__ == "__main__":
    main()
