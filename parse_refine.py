import sys; sys.path.insert(0,'/home/madiyar/uc/scripts')
from local_refine import pocket_residues, pose_contacts, OUTDIR
from pathlib import Path
res = pocket_residues()
rows = []
for p in sorted(Path('/home/madiyar/uc/screening/refine_local').glob('*.pdbqt')):
    if p.name.endswith('.ligand.pdbqt'): continue
    tag = p.stem.split('__')[0]
    vina = [l for l in p.read_text().split('\n') if l.startswith('REMARK VINA RESULT')]
    if not vina: continue
    best = float(vina[0].split()[3])
    allp = [round(float(l.split()[3]), 2) for l in vina]
    c = pose_contacts(p, res)
    rows.append((tag, best, allp, c))
rows.sort(key=lambda r: -r[1])
with open('/home/madiyar/uc/screening/refine_local/refine_results.csv', 'w') as fh:
    fh.write("compound,best_dG,all_pose_dG,contacts\n")
    for tag, best, allp, c in rows:
        fh.write(f"{tag},{best:.2f},\"{allp}\",\"{';'.join(k+':'+str(v) for k, v in list(c.items())[:10])}\"\n")
print("%-16s%8s   nearest pocket contacts (<4.5 A)" % ("compound", "dG"))
for tag, best, allp, c in rows:
    cs = ", ".join(k.split()[0] + ":" + str(v) for k, v in list(c.items())[:6])
    print("%-16s%8.2f   %s" % (tag, best, cs))
print("\nresidues contacted per compound:", [len(r[3]) for r in rows])
