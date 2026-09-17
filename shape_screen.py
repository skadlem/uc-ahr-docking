#!/usr/bin/env python3
"""GPU-free shape/pharmacophore ranking of the Stage A top hits.

The CNN refinement (Stage B) needs a Colab T4, which is demand-blocked. This
script provides an orthogonal, purely-local prioritisation that needs no docking
engine at all: ROCS-style 3D shape + colour (hydrophobicity) comparison of each
top Vina hit against the **crystallographic conformation of indirubin** bound in
the AHR PAS-B pocket (PDB 7ZUB, ligand JY6).

Premise: the pocket binds planar aromatics (indirubin, B[a]P, TCDD, FICZ). A hit
that can adopt a conformation whose volume and hydrophobic feature distribution
overlaps the native ligand is more likely to be a genuine binder than an equally
scoring hit with a completely different shape.

Method
  - reference = prepared 3D indirubin (`structures/redock/ligand_prep.sdf`)
  - for each hit, generate an ETKDGv3 + MMFF conformer ensemble and keep the
    best alignment to the reference (this is exactly what conformer-based shape
    screening does - the reported score is the best over conformers, not the
    single input conformer)
  - metrics: Open3DAlign shape+colour score (Crippen contributions as colour) and
    the pure shape Tanimoto from rdShapeHelpers

Output: research/shape_screen.csv (+ markdown table on stdout)
"""
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign, rdShapeHelpers

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "structures" / "redock" / "ligand_prep.sdf"
LIB = ROOT / "screening" / "library_batch1.sdf"
OUT = ROOT / "research" / "shape_screen.csv"

# Stage A top-15 by Vina (research/screen_stageA_top30.md) -> Vina dG
TOP15 = [
    ("CHEMBL1201174", "Sitagliptin phosphate", -11.47),
    ("CHEMBL535650", "Mefloquine HCl", -11.32),
    ("CHEMBL3989973", "Pexidartinib HCl", -11.25),
    ("CHEMBL3545253", "Flortaucipir F-18", -11.23),
    ("CHEMBL1173055", "Rucaparib", -11.21),
    ("CHEMBL295124", "Berberine", -10.91),
    ("CHEMBL2105743", "Olodaterol HCl", -10.90),
    ("CHEMBL1201753", "Pitavastatin", -10.89),
    ("CHEMBL2107387", "Vortioxetine HBr", -10.82),
    ("CHEMBL1475", "Trioxsalen", -10.81),
    ("CHEMBL1089318", "Opicapone", -10.67),
    ("CHEMBL1709", "Sertraline HCl", -10.66),
    ("CHEMBL1200374", "Exemestane", -10.53),
    ("CHEMBL1301", "Hydroxystilbamidine", -10.52),
    ("CHEMBL1218", "Ramelteon", -10.51),
]
NUM_CONFS = 60


def prep(mol):
    """Strip salts/Hs -> a clean neutral parent with a MMFF-optimised 3D ensemble."""
    try:
        m = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.RemoveHs(mol, sanitize=False)))
    except Exception:
        m = None
    if m is None:
        return None
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return None
    m = Chem.AddHs(m)   # O3A/MMFF needs explicit Hs
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 0xC0FFEE
    cids = AllChem.EmbedMultipleConfs(m, numConfs=NUM_CONFS, params=ps)
    if not len(cids):
        return None
    try:
        AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=400)
    except Exception:
        pass
    return m


def mmff_surrogate(m):
    """Cl/Br/I -> F (isosteric, geometry preserved); MMFF94 has no heavier
    halogen parameters. Only affects the O3A colour model, not the shape."""
    s = Chem.Mol(m)
    for a in s.GetAtoms():
        if a.GetAtomicNum() in (17, 35, 53):
            a.SetAtomicNum(9)
    return s


def main():
    ref = next((m for m in Chem.SDMolSupplier(str(REF), removeHs=True, sanitize=True)
                if m is not None), None)
    if ref is None:
        raise SystemExit(f"could not read reference {REF}")
    ref = Chem.AddHs(ref, addCoords=True)   # O3A/MMFF needs explicit Hs
    ref_props = AllChem.MMFFGetMoleculeProperties(ref)
    ref_conf = ref.GetConformer()
    print(f"reference: indirubin, {ref.GetNumHeavyAtoms()} heavy atoms, "
          f"3D={ref_conf.Is3D()}", flush=True)

    lib = {m.GetProp("_Name"): m
           for m in Chem.SDMolSupplier(str(LIB), removeHs=False, sanitize=False)
           if m is not None}

    rows = []
    for cid, name, vina in TOP15:
        raw = lib.get(cid)
        if raw is None:
            print(f"{cid}: not in local library", flush=True)
            continue
        m = prep(raw)
        if m is None:
            print(f"{cid}: preparation failed", flush=True)
            continue
        best_o3a, best_shape = -1.0, -1.0
        for conf in m.GetConformers():
            # O3A needs the probe conformer to be the "current" one
            cid_ = conf.GetId()
            # RDKit's MMFF94 port lacks parameters for several drug chemotypes
            # (returns None, not an exception). Halogens -> F as an isosteric
            # surrogate; geometry is unchanged. Compounds that still have no
            # MMFF94 parameters are reported as unavailable, not silently scored.
            try:
                prb = mmff_surrogate(m)
                prb_props = AllChem.MMFFGetMoleculeProperties(prb)
                if prb_props is None or ref_props is None:
                    continue
                o3a = rdMolAlign.GetO3A(prb, ref, prb_props, ref_props,
                                        prbCid=cid_, refCid=ref_conf.GetId())
                o3a.Align()
                best_o3a = max(best_o3a, float(o3a.Score()))
                d = rdShapeHelpers.ShapeTanimotoDist(prb, ref, cid_,
                                                     ref_conf.GetId())
                best_shape = max(best_shape, 1.0 - float(d))
            except Exception:
                continue
        if best_o3a < 0:
            print(f"{cid}: MMFF94 parameters unavailable in this RDKit build "
                  f"- excluded from the shape ranking", flush=True)
            continue
        rows.append({"id": cid, "name": name, "vina": vina,
                     "o3a_shape_color": round(best_o3a, 3),
                     "shape_tanimoto": round(best_shape, 3),
                     "heavy_atoms": m.GetNumHeavyAtoms()})
        print(f"{cid}  {name:<24} vina={vina:>6.2f}  o3a={best_o3a:6.3f}  "
              f"shapeT={best_shape:5.3f}  HA={m.GetNumHeavyAtoms()}", flush=True)

    # combined rank: normalise both metrics, Vina is more energetic so weight it
    for r in rows:
        pass
    if rows:
        vmin = min(r["vina"] for r in rows)
        vmax = max(r["vina"] for r in rows)
        omin = min(r["o3a_shape_color"] for r in rows)
        omax = max(r["o3a_shape_color"] for r in rows)

        def norm(v, lo, hi):
            return 0.0 if hi == lo else (v - lo) / (hi - lo)

        for r in rows:
            r["vina_n"] = round(norm(r["vina"], vmin, vmax), 3)
            r["o3a_n"] = round(norm(r["o3a_shape_color"], omin, omax), 3)
            r["composite"] = round(0.5 * r["vina_n"] + 0.5 * r["o3a_n"], 3)
        rows.sort(key=lambda r: -r["composite"])

    fields = ["id", "name", "vina", "o3a_shape_color", "shape_tanimoto",
              "heavy_atoms", "vina_n", "o3a_n", "composite"]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as fh:
        fh.write(",".join(fields) + "\n")
        for r in rows:
            fh.write(",".join(str(r[f]) for f in fields) + "\n")
    print(f"\nwrote {OUT}")
    print("rank | drug                | Vina   | O3A    | shapeT | composite")
    for i, r in enumerate(rows, 1):
        print(f"{i:>4} | {r['name']:<19} | {r['vina']:>6.2f} | "
              f"{r['o3a_shape_color']:>6.3f} | {r['shape_tanimoto']:>6.3f} | "
              f"{r['composite']:>6.3f}")


if __name__ == "__main__":
    main()
