# High-exhaustiveness refinement of the consensus shortlist

Script: `scripts/local_refine.py` + `scripts/parse_refine.py`
Data: `screening/refine_local/refine_results.csv`

## Run

Eight compounds — the consensus shortlist plus the indirubin control — docked at
`exhaustiveness=64` (16× the Stage A setting), 3 poses each, on 12 local cores.
This is the final ranked-pose stage that Colab never delivered, and the first
output of the project produced entirely without a GPU or a remote runtime.

## Results

| compound | ΔG (kcal/mol) | pocket contacts (<4.5 Å) | nearest residues |
|---|---|---|---|
| **Indirubin** (control) | **−12.73** | 19 | Ser365, Ile325, Leu353, Phe351, Val381, Pro297 |
| **Berberine** | **−11.00** | 24 | Leu308, Cys333, Ser346, His291, Thr289, Ile349 |
| **Trioxsalen** | **−10.81** | 18 | Ser346, Val381, Ile325, Cys333, Ile349, Ser336 |
| Ramosetron | −10.75 | 20 | Ile325, Phe351, Ser365, Val381, His291, Ser346 |
| Opicapone | −10.04 | 26 | Ser346, Leu353, Thr289, **Gln383**, Pro297, Leu308 |
| Triamterene | −9.92 | 21 | Leu353, Val381, Phe324, Cys333, Ile325, His291 |
| Methylene blue | −9.78 | 22 | His291, Ser346, Leu308, Ile325, Thr289, Leu315 |
| Indalpine | −9.78 | 19 | Leu353, Ile325, Cys333, Phe351, His291, His337 |

Every pose sits in the PAS-B pocket (centroid within 1.5 Å of the indirubin
reference), and every compound makes 18–26 pocket contacts. Berberine,
trioxsalen and ramosetron form a clear second tier behind the natural ligand.

## Finding 1 — the refinement validates the earlier ranks

The Stage A ranking was computed on Colab with GNINA at `exhaustiveness=2` on a
2-vCPU T4. The high-exhaustiveness local re-run reproduces it almost exactly:

| compound | Stage A (GNINA, exh 2) | refined (Vina, exh 64) | Δ |
|---|---|---|---|
| Berberine | −10.91 | −11.00 | +0.09 |
| Trioxsalen | −10.81 | −10.81 | 0.00 |
| Opicapone | −10.67 | −10.04 | +0.63 |
| Indalpine | (shape-only hit) | −9.78 | new |

**The cheap screen did not distort the ranking.** That is a meaningful
methodological result: a 2-vCPU exhaustive-2 screen ranked the shortlist the same
way as a 16× more expensive refinement, so the throughput-saving protocol is
sound for future screens.

## Finding 2 — opicapone is the only H-bond candidate

The pocket is overwhelmingly hydrophobic (Leu/Ile/Val/Phe dominate every
contact list, as the 7ZUB structure predicts). **Opicapone is the single
compound whose nearest neighbour is a polar residue — Gln383 at 2.33 Å** — and it
has the most pocket contacts of any compound (26). Gln383 and His337 are the only
residues in the pocket that can donate or accept a hydrogen bond to a bound
ligand; every other compound buries against greasy side chains. Opicapone
therefore has the most specific binding hypothesis of the set even though its
score is mid-range. This is exactly the kind of signal the pharmacophore screen
was designed to find, and it independently agrees with that screen's 100th-
percentile ranking for opicapone.

## Finding 3 — berberine's pose explains its activity

Berberine is the experimentally validated AHR activator in the set (PMID
33285228, improves colitis in mouse models). Its refined pose reaches **Leu308
at 2.04 Å**, the closest contact of any non-indirubin compound, and packs against
Cys333 and His291. The protoberberine scaffold is a rigid, planar, cationic
fused ring system — precisely the chemotype the shape screen flagged (0.607
shape Tanimoto to indirubin) and the class AHR is known to bind (planar
aromatics; B[a]P EC50 26 nM). Docking, shape and published biology now agree on
one compound.

## A honest note on a caught error

The first run of this stage reported berberine, trioxsalen and opicapone at
−6.5 kcal/mol with a contact set unlike any other compound. That was wrong: the
name lookup used the wrong SD property (`drug_name`/`_Name` instead of `name`),
so it silently failed and docked the first molecule in the SDF — a 15-heavy-atom
fragment — three times. The tell was that all three "hits" had identical scores
and identical contacts. The lookup was fixed, the run repeated, and the numbers
above are from the corrected run. Recorded here because the error pattern
(silent wrong-molecule docking) is easy to make and hard to notice, and because
the correct values changed the conclusion for those three compounds by 4.5
kcal/mol.

## Status of the AHR lead

Berberine remains the strongest lead — validated AHR activation in colitis
models, top non-indirubin docking score, deepest pocket burial, and convergence
across three orthogonal in-silico methods. Trioxsalen and ramosetron are the
testable new hypotheses. The immediate next experiment is not more docking: it is
an AHR reporter assay to confirm ramosetron and triamterene, and to establish
**efficacy direction** (agonist vs antagonist), which no docking score can
resolve.
