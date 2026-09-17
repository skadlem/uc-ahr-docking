# Local AutoDock-Vina redocking gate — PASS

## The unlock

GNINA cannot run on this box (every release binary dynamically links
`libcudart.so.12`, `libcusparse.so.12`, `libcublas.so.12`), and the free Colab T4
was demand-blocked for the entire session, stalling all docking work. That
bottleneck is now removed: **AutoDock-Vina's official pip wheel (`vina` 1.2.7,
pure CPU, no CUDA)** installs and runs locally on 12 cores, so the entire screen
can be executed on this machine — no Colab, no GPU, no queue.

```
uv pip install --python .venv/bin/python vina   # vina 1.2.7, 6.5 MB, no CUDA
```

## Gate protocol

Same logic as the GNINA gate: dock crystallographic indirubinin (7ZUB chain D,
ligand JY6, 20 heavy atoms) back into the AHR PAS-B pocket and measure how
closely the top pose recovers the crystal binding mode. Receptor preparation
mirrors the Colab run: PDB → pH 7.4 protonation (obabel) → rigid PDBQT
(`-xr`) with Gasteiger charges. Box centred on the reference ligand
(160.76, 164.21, 159.80), 22 Å cubic, `exhaustiveness=32`, 10 modes.

Three non-obvious problems had to be solved to get a trustworthy number:

1. **Receptor PDBQT must be rigid.** A plain `obabel -o pdbqt` writes the protein
   as a flexible *ligand* (724 detected torsions, `ROOT`/`BRANCH` tags), which
   Vina rejects: *"Unknown or inappropriate tag found in rigid receptor"*.
   `-xr` disables torsion detection and produces the rigid receptor Vina needs.
2. **Charges and hydrogens are not added by default.** A raw `-xr` conversion
   leaves every charge at 0.000 and drops all hydrogens, which would make the
   force field meaningless. Protonate at pH 7.4 *first*, then convert.
3. **RDKit cannot read Vina's PDBQT output, and obabel round-trips of it destroy
   aromaticity** (benzene rings come back as `[c][c][c]...`, so no substructure
   match to the reference exists). Fix: PDBQT columns 1–54 *are* a legal PDB
   ATOM record — only the trailing vdW/charge/type columns differ. Truncating
   each line to 54 characters makes it parseable by `MolFromPDBBlock` with
   `proximityBonding=True`, after which `AssignBondOrdersFromTemplate` restores
   correct bond orders and GetBestRMS / coordinate comparison work normally.

## Result

| metric | value | criterion |
|---|---|---|
| Top-pose affinity | **−12.59 kcal/mol** | — |
| **Receptor-frame RMSD (no superposition, identity correspondence)** | **0.540 Å** | < 2.0 Å |
| Rigid-translation component | 0.51 Å | — |
| **Internal-geometry RMSD (translation removed)** | **0.181 Å** | — |
| Per-atom max deviation | 0.671 Å | — |

Every one of the 20 heavy atoms lands within 0.67 Å of its crystallographic
position, and the residual is almost entirely a single 0.51 Å rigid translation
of an otherwise perfectly reproduced geometry (0.18 Å internal RMSD). The docked
pose is the crystal pose.

**GATE: PASS.** For reference the GNINA/Colab gate returned 0.49 Å; the local
Vina receptor-frame result of 0.54 Å is equivalent within the ~0.5 Å resolution
of the 2.9 Å cryo-EM map.

### Validation of the measurement itself

RMSD against a symmetric ligand is easy to fool, so the metric was validated with
negative controls before being trusted: self-vs-self 0.000 Å; Gaussian noise
σ=1.5 Å → 2.678 Å (correctly large). An early hand-rolled Hungarian matcher gave
a *false* 0.106 Å for a 4 Å displacement because it optimised the atom
correspondence on an already-superposed frame; that approach was discarded in
favour of RDKit's `GetBestRMS` plus an explicit receptor-frame identity check, and
the per-atom displacement analysis above is what the conclusion rests on.

Consequence: **the 20-compound Stage D validation library, the missing chunks 4–5
of Stage A, and any refinement can all be run locally now.**

Reproduce: `scripts/local_vina_redock.py` (docking) + `scripts/docking_rmsd_local.py` (RMSD).
