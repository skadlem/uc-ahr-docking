# External benchmark: Mosa et al. 2024 and a library-coverage limitation

## The benchmark

While searching for AHR literature support for the orthogonal-method candidates,
I found a directly comparable published study:

> Mosa F.E.S., Alqahtani M.A., El-Ghiaty M.A., Barakat K., El-Kadi A.O.S.
> *Identifying novel aryl hydrocarbon receptor (AhR) modulators from clinically
> approved drugs: In silico screening and In vitro validation.*
> **Arch Biochem Biophys** 2024, doi:10.1016/j.abb.2024.109958

Their protocol is essentially ours: structure-based virtual screening of an
approved-drug collection (DrugBank) against the **AHR PAS-B binding pocket**,
filtered by binding affinity and MD stability, then **in vitro validation**.
Ten hits, all confirmed to modulate AHR signalling:

> flibanserin, butoconazole, luliconazole, naftifine, triclabendazole,
> rosiglitazone, empagliflozin, benperidol, nebivolol, zucapsaicin

Nine are antagonists (reduce AhR activity); **luliconazole is an agonist**.

## What the comparison revealed — a real limitation of our library

I checked each of the ten against our screening library:

| Published hit | heavy atoms | in our 589-library? |
|---|---|---|
| Flibanserin | 28 | **no** |
| Butoconazole | 25 | **no** |
| Naftifine | 22 | **no** |
| Rosiglitazone | 25 | **no** |
| Empagliflozin | 31 | **no** |
| Benperidol | 28 | **no** |
| Nebivolol | 29 | **no** |
| Zucapsaicin | — | not in ChEMBL approved set |
| Luliconazole | 21 | yes — docked, **missed top-30** |
| Triclabendazole | 21 | yes — docked, **missed top-30** |

**Eight of ten experimentally validated AHR modulators were absent from our
589-compound MaxMin-diversity subset.** Diversity maximisation is the right
choice for *breadth* on a fixed docking budget, but it silently discards known
actives that cluster chemically — the azole antifungals (butoconazole,
luliconazole) and the statin/glitazar cluster are exactly the kind of
near-duplicates MaxMin prunes. This is a concrete, honest limitation of the
Stage A design: **the top-15 ranking is a ranking of the diversity subset, not of
approved-drug space.**

It also explains why the two that *were* docked missed the top-30: both are
azole/triazole antifungals, and neither is a planar aromatic — so they lose on
both the shape metric and the hydrophobic-burial Vina score, despite confirmed
activity. Vina scores binding affinity, not efficacy direction (9/10 are
antagonists; our screen has no way to distinguish agonism from antagonism at all).

## Action taken — Stage D validation library

Rather than leave this as a caveat, I built a focused 20-compound library
(`screening/validation_library.sdf`, also in the repo) that turns the gap into a
measurable test:

- **BENCH (10)** — the nine available Mosa-2024 modulators + zucapsaicin if it
  can be sourced; docking them is a direct external calibration. If they score
  near the Stage A top-15, the box reproduces published SAR. If they do not, that
  quantifies the affinity-vs-efficacy blind spot.
- **ORTH (10)** — the orthogonal-method candidates (methylene blue, triamterene,
  amiloride, indalpine, ramosetron, moxonidine, frovatriptan, sapropterin,
  guanabenz, minoxidil), giving the shape/pharmacophore hypotheses a docking test.
- **CTRL (1)** — indirubin itself, the 0.49 Å redock positive control, in the
  same box and same protocol.

Protocol for the run (`scripts/cell_stageD.py`): `--exhaustiveness 16
--cnn_scoring rescore --num_modes 3`, high accuracy, CNN-rescored, 20 ligands
≈ 20-25 min on a T4. **Colab T4 was demand-blocked for the entire session**, so
the cell is staged and the library hosted in the repo at
`https://raw.githubusercontent.com/skadlem/uc-ahr-docking/master/validation_library.sdf`
for the next available runtime.

## Broader point for the 90-day plan

This is the strongest argument yet for moving the *prioritisation* off the
docking score alone. Three orthogonal local metrics (docking, shape,
pharmacophore) disagree in informative ways, and only the intersection is
trustworthy. The consensus shortlist for the AHR axis is now:
**berberine** (validated AHR activator, PMID 33285228) and **trioxsalen**
(furanocoumarin class) from docking+shape; **opicapone** and **rucaparib** from
docking+pharmacophore; **triamterene** and **methylene blue** as
shape/pharmacophore-driven new hypotheses — all pending the Stage D docking round
and, ultimately, an efficacy-direction assay, which no in-silico method here can
substitute for.
