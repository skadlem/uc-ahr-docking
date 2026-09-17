# 5. STAGE B - CNN refinement of the top Vina hits
# Two environment quirks, both worked around here:
#  - with --cnn_scoring none GNINA writes the Vina score as <minimizedAffinity>,
#    not <Affinity>;
#  - RDKit's SDMolSupplier mis-reads GNINA's .sdf.gz in this Colab image (it
#    sees the whole file as one unparseable record), so every file is
#    decompressed to plain .sdf first. A checkpoint truncated by an interrupted
#    run fails at gzip-open and is skipped instead of aborting the stage.
from rdkit import Chem
import gzip, glob, os, subprocess, time, pandas as pd

def vina(m):
    for k in ('minimizedAffinity', 'Affinity'):
        if m.HasProp(k):
            try:
                return float(m.GetProp(k))
            except ValueError:
                pass
    return None

def read_stageA(f):
    """Decompress + parse one checkpoint. A checkpoint truncated by an
    interrupted run unzips to a short SDF; we keep whatever parsed before the
    truncation point instead of letting the supplier exception kill the stage."""
    out = f.replace('.sdf.gz', '.sdf')
    try:
        with gzip.open(f, 'rt', errors='replace') as fh:
            open(out, 'w').write(fh.read())
    except Exception as e:
        print(f'{f}: corrupt gzip, skipped ({e})', flush=True)
        return []
    mols = []
    try:
        for m in Chem.SDMolSupplier(out, removeHs=True, sanitize=False):
            if m is not None:
                mols.append(m)
    except Exception as e:
        print(f'{f}: truncated SDF, kept {len(mols)} usable poses', flush=True)
    return mols

best = {}
for f in sorted(glob.glob('stageA_*.sdf.gz')):
    mols = read_stageA(f)
    if not mols:
        continue
    n = 0
    for m in mols:
        if m is None:
            continue
        aff = vina(m)
        if aff is None:
            continue
        cid = m.GetProp('_Name') or (m.GetProp('chembl_id') if m.HasProp('chembl_id') else '')
        nm = m.GetProp('name') if m.HasProp('name') else cid
        if not cid:
            continue
        n += 1
        if cid not in best or aff < best[cid][0]:
            best[cid] = (aff, nm)
    print(f'{f}: {n} scored', flush=True)

print(f'{len(best)} unique compounds with a Vina score in stage A', flush=True)
pd.DataFrame([(cid, aff, nm) for cid, (aff, nm) in best.items()],
             columns=['id', 'vina_affinity', 'name']).sort_values(
             'vina_affinity').to_csv('stageA_results.csv', index=False)
print('wrote stageA_results.csv', flush=True)

lib = {}
for m in Chem.SDMolSupplier('library.sdf', removeHs=True, sanitize=False):
    if m is not None:
        lib[m.GetProp('_Name')] = m
print(f'library index: {len(lib)} structures', flush=True)

# If the stage A checkpoints are absent (e.g. a fresh runtime after a
# disconnect), fall back to the top-15 Vina hits recorded in
# research/screen_stageA_top30.md and go straight to CNN refinement.
KNOWN_TOP15 = ['CHEMBL1201174', 'CHEMBL535650', 'CHEMBL3989973', 'CHEMBL3545253',
               'CHEMBL1173055', 'CHEMBL295124', 'CHEMBL2105743', 'CHEMBL1201753',
               'CHEMBL2107387', 'CHEMBL1475', 'CHEMBL1089318', 'CHEMBL1709',
               'CHEMBL1200374', 'CHEMBL1301', 'CHEMBL1218']

TOPN = 30
ranked = sorted(best.items(), key=lambda kv: kv[1][0])[:TOPN]
if not ranked:
    print('no stage A checkpoints found - using recorded top-15 Vina hits', flush=True)
    ranked = [(cid, (0.0, '')) for cid in KNOWN_TOP15]
w = Chem.SDWriter('shortlist.sdf')
n_ok = 0
for cid, (aff, nm) in ranked:
    if cid in lib:
        w.write(lib[cid])
        n_ok += 1
w.close()
print(f'shortlist: {n_ok}/{len(ranked)} top Vina hits -> shortlist.sdf', flush=True)
for cid, (aff, nm) in ranked[:15]:
    print(f'  {cid}  {nm[:26]:<26} vina={aff:.2f}', flush=True)

t = time.time()
r = subprocess.run(['./gnina', '-r', 'receptor_prep.pdb', '-l', 'shortlist.sdf',
                    '--autobox_ligand', 'indirubin.pdb', '--autobox_add', '5',
                    '--cnn', 'fast', '--cnn_scoring', 'rescore',
                    '--exhaustiveness', '12', '--num_modes', '3',
                    '-o', 'refined.sdf.gz'],
                   capture_output=True, text=True, timeout=5400)
print(f'STAGE B rc={r.returncode} in {(time.time()-t)/60:.1f} min', flush=True)
if r.returncode != 0:
    print((r.stderr or '')[-500:], flush=True)
