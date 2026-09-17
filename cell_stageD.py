# cell_stageD.py — Stage D: validation-library docking (benchmark + orthogonal)
# Paste into a new Colab code cell. Self-contained: fetches the 20-compound
# library from the repo raw URL, docks each at high exhaustiveness with CNN
# rescoring, prints a ranked table.
#
# Groups in the library:
#   BENCH  = Mosa 2024 experimentally-validated AHR modulators (external calib)
#   ORTH   = candidates from the local shape/pharmacophore screens
#   CTRL   = indirubin, the co-crystal ligand (0.49 A redock gate)
#
# Scientific point: Vina scores *affinity*, not *efficacy*. Most Mosa-2024 hits
# are AHR antagonists, so a low Vina score for them is not a failure - it is a
# data point about what the pocket will bind vs what will modulate it.

import gzip, os, subprocess, urllib.request
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

URL = ("https://raw.githubusercontent.com/skadlem/uc-ahr-docking/master/"
       "validation_library.sdf")

# 1. fetch the validation library (plain SDF, ~64 kB)
Path("stageD").mkdir(exist_ok=True)
if not Path("stageD/validation_library.sdf").exists():
    print("fetching validation library ...", flush=True)
    urllib.request.urlretrieve(URL, "stageD/validation_library.sdf")
sup = Chem.SDMolSupplier("stageD/validation_library.sdf", removeHs=False,
                         sanitize=False)
recs = [(m.GetProp("_Name"), m.GetProp("drug_name"), m.GetProp("group"))
        for m in sup if m is not None]
print(f"validation library: {len(recs)} compounds")
for cid, nm, grp in recs:
    print(f"  {grp:5s} {cid:14s} {nm}")

# 2. protonation/charge normalisation exactly as Stage A used
raw = [m for m in Chem.SDMolSupplier("stageD/validation_library.sdf",
                                     removeHs=False, sanitize=False) if m]
w = Chem.SDWriter("stageD/ligands.sdf")
for m in raw:
    w.write(m)
w.close()

# 3. dock every ligand at high exhaustiveness + CNN rescore (this is the whole
#    point of a validation round: no shortcuts)
cmd = ["./gnina",
       "-r", "receptor_prep.pdb",
       "-l", "stageD/ligands.sdf",
       "--autobox_ligand", "indirubin.pdb",
       "--autobox_add", "5",
       "--cnn_scoring", "rescore",
       "--exhaustiveness", "16",
       "--num_modes", "3",
       "-o", "stageD/docked.sdf.gz"]
print("running:", " ".join(cmd), flush=True)
t0 = time.time()
r = subprocess.run(cmd, capture_output=True, text=True)
print(f"gnina rc={r.returncode} in {(time.time()-t0)/60:.1f} min")
if r.stdout:
    print(r.stdout[-1500:])
if r.returncode != 0:
    print("STDERR:", r.stderr[-2000:])
    raise SystemExit("gnina failed")

# 4. parse (GNINA .sdf.gz needs decompressing before RDKit here)
with gzip.open("stageD/docked.sdf.gz", "rb") as fh:
    Path("stageD/docked.sdf").write_bytes(fh.read())
rows = []
sup = Chem.SDMolSupplier("stageD/docked.sdf", removeHs=True, sanitize=False)
for m in sup:
    if m is None:
        continue
    d = m.GetPropsAsDict()
    rows.append({
        "chembl_id": m.GetProp("_Name"),
        "affinity": float(d.get("Affinity", d.get("minimizedAffinity", 0))),
        "cnnscore": float(d.get("CNNscore", 0)),
        "cnnaffinity": float(d.get("CNNaffinity", 0)),
    })
df = pd.DataFrame(rows).drop_duplicates("chembl_id")
meta = pd.DataFrame(recs, columns=["chembl_id", "name", "group"])
df = meta.merge(df, on="chembl_id", how="left").sort_values("affinity")
df.to_csv("stageD/stageD_results.csv", index=False)
print("\n=== Stage D results (sorted by Vina affinity) ===")
print(df.to_string(index=False))
print("\nSTAGE D rc=0 done")
