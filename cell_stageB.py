from rdkit import Chem
import glob, subprocess, time, pandas as pd

best = {}
for f in sorted(glob.glob('stageA_*.sdf.gz')):
    for m in Chem.SDMolSupplier(f, removeHs=True, sanitize=False):
        if m is None:
            continue
        cid = m.GetProp('_Name')
        if not m.HasProp('Affinity'):
            continue
        aff = float(m.GetProp('Affinity'))
        if cid not in best or aff < best[cid][0]:
            best[cid] = (aff, m)

print(f'{len(best)} unique compounds scored in stage A')
rows = [(cid, aff, m.GetProp('name') if m.HasProp('name') else cid) for cid, (aff, m) in best.items()]
pd.DataFrame(rows, columns=['id', 'vina_affinity', 'name']).sort_values('vina_affinity').to_csv('stageA_results.csv', index=False)
print('wrote stageA_results.csv (all compounds, Vina affinity)')

TOPN = 30
ranked = sorted(best.items(), key=lambda kv: kv[1][0])[:TOPN]
w = Chem.SDWriter('shortlist.sdf')
for cid, (aff, m) in ranked:
    w.write(m)
w.close()
print(f'shortlist: {len(ranked)} top Vina hits -> shortlist.sdf', flush=True)

t = time.time()
r = subprocess.run(['./gnina', '-r', 'receptor_prep.pdb', '-l', 'shortlist.sdf',
                    '--autobox_ligand', 'indirubin.pdb', '--autobox_add', '5',
                    '--cnn', 'fast', '--cnn_scoring', 'rescore',
                    '--exhaustiveness', '12', '--num_modes', '3',
                    '-o', 'refined.sdf.gz'],
                   capture_output=True, text=True, timeout=5400)
print(f'STAGE B rc={r.returncode} in {(time.time()-t)/60:.1f} min')
if r.returncode != 0:
    print((r.stderr or '')[-400:])
