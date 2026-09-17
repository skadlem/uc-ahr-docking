#!/usr/bin/env python3
"""Docking-independent virtual screen of the FULL Stage A library.

Motivation
  Three independent lines of evidence now say the AHR PAS-B pocket selects
  *planar aromatic* ligands:
    1. the co-crystal ligand indirubin, B[a]P, TCDD, FICZ are all flat aromatics;
    2. the Stage A docking top hits separate into planar aromatics and flexible
       non-planar high-scoring decoys;
    3. the O3A shape screen (`shape_screen.py`) ranked the top-15 by shape
       Tanimoto and the four most indirubin-like are all planar aromatics.

  This turns that observation into a *screening criterion*: score every compound
  in the 424-member Stage A library (<=34 heavy atoms) by how well it can match
  indirubin's bound shape, and by how planar it is - both computed locally with
  RDKit and no docking engine. The intersection with the Vina ranking is the
  consensus shortlist.

Metrics per compound
  - shape_tanimoto : best over a 25-conformer ETKDGv3+MMFF ensemble vs the
                     crystallographic indirubin pose (Open3DAlign)
  - planarity      : RMSD (A) of heavy atoms from their best-fit plane
  - aromatic_rings : RDKit SSSR aromatic-ring count
  - fused_rings    : count of rings sharing >=2 bonds (ring-system size)

Resumable: appends one CSV row per compound, so an interrupted run just resumes.
Parallel over all cores.
"""
import csv
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign, rdShapeHelpers
from rdkit.Geometry import Point3D

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
LIB = ROOT / "screening" / "library_batch1.sdf"
REF = ROOT / "structures" / "redock" / "ligand_prep.sdf"
OUT = ROOT / "screening" / "planarity_screen.csv"
MAX_HA = 34        # matches the Stage A filter
NUM_CONFS = 25

_REF = None


def _load_ref():
    global _REF
    if _REF is None:
        m = next((x for x in Chem.SDMolSupplier(str(REF)) if x is not None), None)
        _REF = Chem.AddHs(m, addCoords=True)
    return _REF


def planarity(mol):
    """RMSD of heavy atoms from their least-squares plane (smaller = flatter)."""
    conf = mol.GetConformer()
    pts = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    n = len(pts)
    if n < 3:
        return 0.0
    # least-squares plane normal via SVD on centred coordinates
    import numpy as np
    A = np.array([[p.x, p.y, p.z] for p in pts])
    ctr = A.mean(axis=0)
    _, _, Vt = np.linalg.svd(A - ctr)
    normal = Vt[2]
    dist = (A - ctr) @ normal
    return float(float((dist ** 2).sum() / n) ** 0.5)


def fused_ring_count(mol):
    rings = [set(r) for r in Chem.GetSymmSSSR(mol)]
    fused = 0
    for i in range(len(rings)):
        for j in range(i + 1, len(rings)):
            if len(rings[i] & rings[j]) >= 2:
                fused += 1
    return fused


def mmff_surrogate(m):
    """Cl/Br/I -> F (isosteric, geometry preserved) - MMFF94 has no heavier
    halogen parameters. Only affects the O3A colour model, never the shape."""
    s = Chem.Mol(m)
    for a in s.GetAtoms():
        if a.GetAtomicNum() in (17, 35, 53):
            a.SetAtomicNum(9)
    return s


def _inertia_shape(prb, ref, prbCid, refCid=0):
    """Shape Tanimoto after a force-field-free alignment.

    For the compounds where RDKit's MMFF94 has no parameters (so no O3A), align
    the probe onto the reference by matching principal axes of inertia and the
    centroid, trying the degenerate sign permutations and keeping the best
    overlap. Deterministic, needs no force field."""
    try:
        def axes(m, cid):
            pts = np.array(m.GetConformer(cid).GetPositions())
            ctr = pts.mean(axis=0)
            C = pts - ctr
            I = np.eye(3) * np.einsum('ij,ij->', C, C) - C.T @ C
            _, vecs = np.linalg.eigh(I)      # columns: principal axes
            return ctr, vecs

        pc, pv = axes(prb, prbCid)
        rc, rv = axes(ref, refCid)
        orig = np.array(prb.GetConformer(prbCid).GetPositions())
        best = -1.0
        for signs in ((1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1),
                      (1, -1, -1), (-1, 1, -1), (-1, -1, 1), (-1, -1, -1)):
            B = pv * np.array(signs)
            if np.linalg.det(B) < 0:
                continue                     # keep right-handed frames only
            R = rv @ B.T
            moved = (orig - pc) @ R.T + rc
            conf = prb.GetConformer(prbCid)
            for k, atom in enumerate(prb.GetAtoms()):
                conf.SetAtomPosition(k, Point3D(*moved[k]))
            try:
                d = rdShapeHelpers.ShapeTanimotoDist(prb, ref, prbCid, refCid)
                best = max(best, 1.0 - float(d))
            except Exception:
                pass
        # restore so the caller's other conformers are unaffected
        conf = prb.GetConformer(prbCid)
        for k, atom in enumerate(prb.GetAtoms()):
            conf.SetAtomPosition(k, Point3D(*orig[k]))
        return best
    except Exception:
        return -1.0


def score_one(args):
    cid, smi = args
    ref = _load_ref()
    try:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        Chem.SanitizeMol(m)
        if m.GetNumHeavyAtoms() > MAX_HA:
            return None
        m = Chem.AddHs(m)
        ps = AllChem.ETKDGv3()
        ps.randomSeed = 1234
        cids = AllChem.EmbedMultipleConfs(m, numConfs=NUM_CONFS, params=ps)
        if not len(cids):
            return None
        try:
            AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=300)
        except Exception:
            try:
                AllChem.UFFOptimizeMoleculeConfs(m, maxIters=300)
            except Exception:
                pass
        # ShapeTanimotoDist does NOT translate the probe itself: the two
        # molecules must already be superimposed, so every conformer is aligned
        # onto indirubin before its shape Tanimoto is read.
        best_shape, best_o3a = -1.0, None
        prb = mmff_surrogate(m)
        prb_props = AllChem.MMFFGetMoleculeProperties(prb)
        ref_props = AllChem.MMFFGetMoleculeProperties(ref)
        for conf in m.GetConformers():
            i = conf.GetId()
            # (a) O3A alignment (MMFF colour model) - also gives the O3A score
            if prb_props is not None and ref_props is not None:
                try:
                    o3a = rdMolAlign.GetO3A(prb, ref, prb_props, ref_props,
                                           prbCid=i, refCid=0)
                    o3a.Align()
                    s_ = float(o3a.Score())
                    best_o3a = s_ if best_o3a is None else max(best_o3a, s_)
                    d = rdShapeHelpers.ShapeTanimotoDist(prb, ref, i, 0)
                    best_shape = max(best_shape, 1.0 - float(d))
                    continue
                except Exception:
                    pass
            # (b) force-field-free fallback: principal-axes + centroid alignment
            best_shape = max(best_shape, _inertia_shape(prb, ref, i, 0))
        heavy = Chem.RemoveHs(m)
        ri = heavy.GetRingInfo()
        arom = sum(1 for r in ri.AtomRings()
                   if all(heavy.GetAtomWithIdx(a).GetIsAromatic() for a in r))
        return {
            "chembl_id": cid,
            "shape_tanimoto": round(best_shape, 4),
            "o3a_score": "" if best_o3a is None else round(best_o3a, 2),
            "planarity": round(planarity(heavy), 3),
            "heavy_atoms": heavy.GetNumHeavyAtoms(),
            "aromatic_rings": arom,
            "fused_rings": fused_ring_count(heavy),
        }
    except Exception:
        return None


def main():
    sup = Chem.SDMolSupplier(str(LIB), removeHs=False, sanitize=False)
    todos = []
    for mol in sup:
        if mol is None:
            continue
        cid = mol.GetProp("_Name")
        try:
            smi = Chem.MolToSmiles(Chem.RemoveHs(mol, sanitize=False))
        except Exception:
            continue
        if Chem.MolFromSmiles(smi).GetNumHeavyAtoms() > MAX_HA:
            continue
        todos.append((cid, smi))

    done = set()
    if OUT.exists():
        with OUT.open() as fh:
            for row in csv.DictReader(fh):
                if row.get("chembl_id"):
                    done.add(row["chembl_id"])
    todo = [t for t in todos if t[0] not in done]
    print(f"library<= {MAX_HA} HA: {len(todos)} | already scored: {len(done)} | "
          f"to do: {len(todo)}", flush=True)

    fields = ["chembl_id", "shape_tanimoto", "o3a_score", "planarity",
              "heavy_atoms", "aromatic_rings", "fused_rings"]
    newfile = not OUT.exists()
    workers = min(12, os.cpu_count() or 4)
    print(f"workers: {workers}", flush=True)
    n_ok, n_fail = 0, 0
    with OUT.open("a", newline="") as fh, ProcessPoolExecutor(max_workers=workers) as ex:
        w = csv.DictWriter(fh, fieldnames=fields)
        if newfile:
            w.writeheader()
        for res in ex.map(score_one, todo, chunksize=2):
            if res is None:
                n_fail += 1
                continue
            w.writerow(res)
            fh.flush()
            n_ok += 1
            if n_ok % 25 == 0:
                print(f"  {n_ok} scored ({n_fail} failed)", flush=True)
    print(f"DONE: {n_ok} scored, {n_fail} failed -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
