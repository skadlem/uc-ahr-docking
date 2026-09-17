import subprocess, time, os
from concurrent.futures import ThreadPoolExecutor
from rdkit import Chem

# Filter to <= 34 heavy atoms: the AHR PAS-B pocket is compact (indirubin = 20 HA),
# so larger molecules rarely fit and dominate runtime on this 2-vCPU instance.
sup = Chem.SDMolSupplier('library.sdf', removeHs=True, sanitize=False)
keep = [m for m in sup if m is not None and m.GetNumHeavyAtoms() <= 34]
w = Chem.SDWriter('library_filtered.sdf')
for m in keep:
    w.write(m)
w.close()
print(f'{len(keep)}/589 pass <=34 heavy atoms', flush=True)

CH, NPROC = 75, 2
nchunks = (len(keep) + CH - 1) // CH
for k in range(nchunks):
    w = Chem.SDWriter(f'lib_{k:02d}.sdf')
    for m in keep[k * CH:(k + 1) * CH]:
        w.write(m)
    w.close()
print(f'{len(keep)} -> {nchunks} chunks of ~{CH}', flush=True)

def run(k):
    out = f'stageA_{k:02d}.sdf.gz'
    env = {**os.environ, 'OMP_NUM_THREADS': '1'}
    t = time.time()
    r = subprocess.run(['./gnina', '-r', 'receptor_prep.pdb', '-l', f'lib_{k:02d}.sdf',
                        '--autobox_ligand', 'indirubin.pdb', '--autobox_add', '5',
                        '--cnn_scoring', 'none', '--exhaustiveness', '2', '--num_modes', '1',
                        '--cpu', '1', '-o', out],
                       capture_output=True, text=True, timeout=5400, env=env)
    return k, r.returncode, time.time() - t, out, (r.stderr or '')[-200:]

t0 = time.time(); ok = 0
with ThreadPoolExecutor(max_workers=NPROC) as ex:
    for k, rc, dt, out, err in ex.map(run, range(nchunks)):
        if rc != 0:
            print(f'chunk {k}: FAILED rc={rc}\n{err}', flush=True)
        else:
            ok += 1
            el = (time.time() - t0) / 60
            print(f'chunk {k}: ok {dt/60:.1f}min | {ok}/{nchunks} | elapsed {el:.1f}min | ETA {el/ok*(nchunks-ok):.1f}min', flush=True)
print(f'STAGE A done: {ok}/{nchunks} chunks in {(time.time()-t0)/60:.1f} min', flush=True)
