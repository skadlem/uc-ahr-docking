#!/usr/bin/env python3
"""Pharmacophore-feature screen of the Stage A library (GPU-free).

Complement to `planarity_screen.py`. Shape similarity asks "does the molecule
fill the same space as indirubin?"; this asks "does it present the same
*functional features* in the same spatial arrangement?" - hydrogen-bond donors
and acceptors, hydrophobic centres, aromatic rings, cations and anions. Two
molecules can share a shape and differ completely in pharmacophore, so a
credible prioritisation needs both.

Method: Gobbi 2D pharmacophore fingerprints (Gobbi & Lee, JCIM 2003,
doi:10.1021/ci0256062) as implemented in RDKit. Each compound is encoded as a
bitstring over all feature triples with their pairwise distances binned; the
similarity to indirubin is the Tanimoto coefficient. 2D (topological) distances
are used, so this is conformation-independent and complements the 3D shape
metric rather than duplicating it.

Output: screening/pharmacophore_screen.csv
"""
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Pharm2D
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
from rdkit.DataStructs import TanimotoSimilarity

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parent.parent
LIB = ROOT / "screening" / "library_batch1.sdf"
REF = ROOT / "structures" / "redock" / "ligand_prep.sdf"
SHAPE = ROOT / "screening" / "planarity_screen_named.csv"
OUT = ROOT / "screening" / "pharmacophore_screen.csv"
MAX_HA = 34

FACTORY = Gobbi_Pharm2D.factory
# distances in topological-bond space, coarse bins as Gobbi published
FACTORY.SetBins([(0, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8),
                 (8, 100)])


def canon_smiles(mol):
    try:
        return Chem.MolToSmiles(Chem.RemoveHs(mol, sanitize=False))
    except Exception:
        return None


def main():
    ref = next((m for m in Chem.SDMolSupplier(str(REF)) if m is not None), None)
    if ref is None:
        raise SystemExit(f"cannot read reference {REF}")
    ref_smi = Chem.MolToSmiles(ref)
    ref_canon = Chem.MolFromSmiles(ref_smi)
    ref_fp = Generate.Gen2DFingerprint(ref_canon, FACTORY)
    ref_bits = set(ref_fp.GetOnBits())
    n_ref = len(ref_bits)
    print(f"reference indirubin: {ref_canon.GetNumHeavyAtoms()} HA, "
          f"{n_ref} pharmacophore bits on", flush=True)

    rows = []
    sup = Chem.SDMolSupplier(str(LIB), removeHs=False, sanitize=False)
    for mol in sup:
        if mol is None:
            continue
        cid = mol.GetProp("_Name")
        smi = canon_smiles(mol)
        if smi is None:
            continue
        try:
            m = Chem.MolFromSmiles(smi)
            if m is None or m.GetNumHeavyAtoms() > MAX_HA:
                continue
            fp = Generate.Gen2DFingerprint(m, FACTORY)
        except Exception:
            continue
        bits = set(fp.GetOnBits())
        if not bits:
            continue
        # Tanimoto over the full bit space
        sim = TanimotoSimilarity(fp, ref_fp)
        # enriched view: of indirubin's features, how many does this compound
        # present at all (feature recall), and how specific is it to indirubin
        recall = len(bits & ref_bits) / n_ref
        rows.append({"chembl_id": cid, "pharm_tanimoto": round(sim, 4),
                     "pharm_recall": round(recall, 4),
                     "n_bits": len(bits)})
        if len(rows) % 100 == 0:
            print(f"  {len(rows)} done", flush=True)

    df = pd.DataFrame(rows)
    # merge the shape/planarity screen for a joint prioritisation
    if SHAPE.exists():
        shape = pd.read_csv(SHAPE, na_values=[""])
        df = df.merge(shape, on="chembl_id", how="inner")
        # rank-normalise each metric so they are comparable, then combine
        for col in ("shape_tanimoto", "pharm_tanimoto", "pharm_recall"):
            df[col + "_r"] = df[col].rank(pct=True)
        # shape carries the pocket-volume evidence; pharmacophore recall checks
        # the H-bond machinery is present; weight shape slightly higher
        df["consensus"] = ((0.45 * df.shape_tanimoto_r
                            + 0.30 * df.pharm_tanimoto_r
                            + 0.25 * df.pharm_recall_r).round(4))
        df = df.sort_values("consensus", ascending=False)

    df.to_csv(OUT, index=False)
    print(f"\nwrote {OUT} ({len(df)} compounds)")
    show = [c for c in ("chembl_id", "name", "consensus", "shape_tanimoto",
                        "pharm_tanimoto", "pharm_recall", "planarity",
                        "aromatic_rings") if c in df.columns]
    print(df.head(25)[show].to_string(index=False))


if __name__ == "__main__":
    main()
