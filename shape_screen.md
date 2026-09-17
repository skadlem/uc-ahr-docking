# Shape/pharmacophore screen of the Stage A top-15 (GPU-free)

Script: `scripts/shape_screen.py` · Data: `research/shape_screen.csv`

## Why this analysis exists

The Stage B CNN refinement needs a Colab T4, which was demand-blocked after the
free-tier runtime disconnected. To keep producing verifiable results on the
local CPU-only box, this is an **orthogonal, docking-independent prioritisation**:
a ROCS-style 3D shape + colour comparison of each top Vina hit against the
**crystallographic conformation of indirubin** as bound in the AHR PAS-B pocket
(PDB 7ZUB, ligand JY6).

Premise, stated as a testable hypothesis: the PAS-B pocket binds *planar
aromatics* (indirubin, benzo[a]pyrene, TCDD, FICZ). A hit that can adopt a
conformation whose volume overlaps the native ligand is more likely to be a
genuine binder than an equally-scored hit of unrelated shape.

## Method

1. **Reference** = prepared 3D indirubin (`structures/redock/ligand_prep.sdf`),
   the exact geometry validated in the 0.49 Å redocking gate.
2. For each hit: strip input 3D, rebuild from SMILES, add explicit Hs, generate a
   **60-conformer ensemble** (ETKDGv3 + MMFF optimisation) — the reported score
   is the *best over conformers*, which is what conformer-based shape screening
   (ROCS) does. A single input conformer would under-report flexible molecules.
3. **Alignment**: RDKit Open3DAlign (`rdMolAlign.GetO3A`), MMFF atom types as the
   "colour" model. `o3a.Align()` superimposes the probe onto indirubin, then the
   score is `o3a.Score()` (shape + colour overlap; higher = better) and the pure
   shape Tanimoto is `1 - rdShapeHelpers.ShapeTanimotoDist`.
4. Halogen → F isosteric surrogate (geometry preserved) because RDKit's MMFF94
   port lacks heavier-halogen parameters. This affects only the colour model.

## Results — shape Tanimoto (the chemically meaningful metric)

Ranked by how well each hit can reproduce indirubin's bound shape:

| Drug | Vina ΔG | shape Tanimoto | heavy atoms |
|---|---|---|---|
| **Trioxsalen** | −10.81 | **0.630** | 17 |
| **Berberine** | −10.91 | **0.607** | 25 |
| Ramelteon | −10.51 | 0.597 | 19 |
| Hydroxystilbamidine | −10.52 | 0.585 | 21 |
| Exemestane | −10.53 | 0.533 | 22 |
| Rucaparib | −11.21 | 0.509 | 24 |
| Flortaucipir F-18 | −11.23 | 0.504 | 20 |
| Sitagliptin | −11.47 | 0.499 | 34 |
| Pitavastatin | −10.89 | 0.486 | 31 |
| Opicapone | −10.67 | 0.458 | 27 |

Reference self-alignment: O3A score 170.32, shape Tanimoto 1.000 (sanity gate).

## Interpretation

**The top four by shape are all planar aromatic systems.** This is an
independent, docking-free corroboration of both the docking results and the
positive-control literature:

- **Trioxsalen** (highest shape, 0.630) is a planar furocoumarin — the same
  chemical class as the classical AHR pro-agonist bergamottin/grapefruit
  furanocoumarins (see `research/positive_controls.md`, class-level support).
- **Berberine** (0.607) is the **recovered positive control**: a published AHR
  activator that improves colitis via microbial tryptophan catabolism
  (PMID 33285228). It now ranks top-2 by shape *and* top-6 by docking.
- Ramelteon (indane) and hydroxystilbamidine (stilbene) are likewise flat,
  rigid, aromatic.

Meanwhile the flexible, non-planar hits (pitavastatin, sitagliptin) score well on
raw Vina but poorly on shape — consistent with a large-scoring-volume but
shape-mismatched fit.

### Caveat on the composite column

`shape_screen.csv` also carries a `composite` = 0.5·normalised(Vina) +
0.5·normalised(O3A score). **The raw O3A score scales with molecular volume**,
so it inflates large molecules (pitavastatin, 31 heavy atoms) relative to the
compact native ligand (20 heavy atoms). The **shape Tanimoto** column is the
size-normalised metric and is the one used for interpretation above; the
composite is reported only for completeness and should not be read as a
prioritisation.

### Compounds excluded

Five of the fifteen could not be scored because RDKit's MMFF94 port returns *no
parameters* for some atom environments even after the halogen surrogate
(`MMFFGetMoleculeProperties` yields `None` rather than raising):

Mefloquine, Pexidartinib, Olodaterol, Vortioxetine, Sertraline.

Olodaterol has no halogen at all, so this is a general MMFF94 coverage gap in
this build (RDKit 2026.03.6), not a halogen issue. These five remain ranked by
their Stage A Vina score only and are flagged in the CSV by absence. Notably
**three of the five are literature-unsupported anyway** (mefloquine, flortaucipir
class, pexidartinib), so the exclusion does not remove any positive control.

## Bottom line

Two independent computational methods (docking, shape) plus the literature now
converge on the same chemotype: **planar aromatics fit the AHR PAS-B pocket**, and
the two best-supported candidates from the approved-drug screen are
**berberine** (docking −10.91, shape 0.607, published AHR activator with colitis
efficacy) and **trioxsalen** (−10.81, 0.630, furanocoumarin class).
