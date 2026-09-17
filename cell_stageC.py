# 7. STAGE C - rank refined hits + pocket contacts
# His337 is the closest residue to the co-bound indirubin (2.8 A);
# Gln383 / Tyr322 / Ser336 are the only plausible H-bond partners.
# NB: CNN scores are near-degenerate across binding modes (see the redocking
# result - the rank-1 pose was a 6.55 A decoy), so we report the BEST pose per
# compound, not the first.
from rdkit import Chem
from Bio.PDB import PDBParser
import numpy as np, pandas as pd

KEY = {'HIS337', 'GLN383', 'TYR322', 'LEU308', 'ILE325', 'LEU353', 'PHE351', 'PHE287'}

rec = PDBParser(QUIET=True).get_structure('R', 'receptor_prep.pdb')[0]
res_atoms = [(res.get_resname() + str(res.id[1]),
              np.array([a.coord for a in res.get_atoms()])) for res in rec.get_residues()]

def g(m, k):
    return float(m.GetProp(k)) if m.HasProp(k) else None

# RDKit mis-reads GNINA's .sdf.gz in this image - decompress first.
import gzip
with gzip.open('refined.sdf.gz', 'rt', errors='replace') as fh:
    open('refined.sdf', 'w').write(fh.read())

best_per = {}
for m in Chem.SDMolSupplier('refined.sdf', removeHs=True, sanitize=False):
    if m is None:
        continue
    cid = m.GetProp('_Name')
    px = m.GetConformer().GetPositions()
    contacts = []
    for name, xyz in res_atoms:
        d = float(np.linalg.norm(xyz - px, axis=1).min())
        if d <= 4.5:
            contacts.append((name, round(d, 1)))
    contacts.sort(key=lambda t: t[1])
    cnn_aff = g(m, 'CNNaffinity')
    info = {
        'id': cid,
        'name': m.GetProp('name') if m.HasProp('name') else cid,
        'cnn_aff': cnn_aff,
        'cnn_pose': g(m, 'CNNscore'),
        'vina': g(m, 'Affinity'),
        'n_contacts': len(contacts),
        'contacts': contacts,
    }
    key = info['cnn_aff'] if info['cnn_aff'] is not None else info['vina']
    if cid not in best_per or key > best_per[cid]['_rankkey']:
        info['_rankkey'] = key
        best_per[cid] = info

ranked = sorted(best_per.values(), key=lambda r: -r['_rankkey'])
print(f'{len(ranked)} refined compounds\n')
hdr = f"{'CHEMBL':<10}{'name':<22}{'CNNaff':>7}{'CNNpose':>8}{'vina':>7}{'#cont':>6}  key-contacts"
print(hdr); print('-' * len(hdr))
for r in ranked[:20]:
    kc = [n for n, d in r['contacts'] if n in KEY]
    ca = f"{r['cnn_aff']:>7.2f}" if r['cnn_aff'] is not None else '     --'
    cp = f"{r['cnn_pose']:>8.3f}" if r['cnn_pose'] is not None else '       --'
    print(f"{r['id']:<10}{r['name'][:21]:<22}{ca}{cp}"
          f"{r['vina']:>7.2f}{r['n_contacts']:>6}  {' '.join(kc[:6])}")

df = pd.DataFrame([{k: v for k, v in r.items() if k not in ('contacts', '_rankkey')}
                   for r in ranked])
df['top_contacts'] = [' '.join(f'{n}({d})' for n, d in r['contacts'][:8]) for r in ranked]
df.to_csv('screen_results.csv', index=False)
print('\nwrote screen_results.csv')
