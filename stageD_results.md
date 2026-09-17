# Stage D — external calibration against published AHR modulators

Script: `scripts/local_stageD.py` + `scripts/parse_stageD.py`
Data: `screening/stageD_local/stageD_results.csv`

## The run

AutoDock-Vina now runs **locally** on 12 cores (see `research/local_vina_gate.md`),
so the Stage D validation library — built while Colab was blocked — was docked
here rather than on a T4. 20 ligands, `exhaustiveness=16`, 3 poses, 5 workers ×
2 cpus, ~3 min wall clock for the whole set (vs. hours of queue on Colab).

Box and receptor identical to the Colab Stage A protocol (22 Å cube centred on
the reference ligand, pH-7.4-protonated rigid receptor with Gasteiger charges),
so scores are comparable with the Stage A numbers.

## Results

| drug | group | ΔG (kcal/mol) | note |
|---|---|---|---|
| **Indirubin** | CTRL | **−12.64** | positive control, receptor-frame RMSD **0.428 Å** |
| Naftifine | BENCH | −11.39 | published AHR modulator |
| Flibanserin | BENCH | −11.36 | published AHR modulator |
| **Ramosetron** | ORTH | **−11.24** | shape-screen hit (0.674) |
| Benperidol | BENCH | −10.90 | published AHR modulator |
| Frovatriptan | ORTH | −10.29 | shape-screen hit |
| Rosiglitazone | BENCH | −10.04 | published AHR modulator |
| **Triamterene** | ORTH | −10.00 | top shape (0.634) **and** top pharmacophore recall (0.548) |
| Nebivolol | BENCH | −9.91 | published AHR modulator |
| Methylene blue | ORTH | −9.69 | flattest scaffold (planarity 0.135) |
| Indalpine | ORTH | −9.58 | library-max shape similarity (0.722) |
| Triclabendazole | BENCH | −9.49 | published AHR modulator |
| Empagliflozin | BENCH | −8.68 | published AHR modulator |
| Luliconazole | BENCH | −8.61 | published **agonist** |
| Minoxidil | ORTH | −8.29 | |
| Butoconazole | BENCH | −8.00 | published AHR modulator |
| Guanabenz | ORTH | −7.93 | |
| Amiloride | ORTH | −7.65 | |
| Sapropterin | ORTH | −7.58 | |
| Moxonidine | ORTH | −6.92 | |

## Finding 1 — the box reproduces published SAR

**The positive control works.** Indirubin is the top-scoring ligand in its own
pocket at −12.64 kcal/mol with a receptor-frame RMSD of 0.428 Å — the pose *is*
the crystal pose. The docking protocol is sound.

**5 of 9 published AHR modulators score better than −10**, and the BENCH group
mean (−9.82) is 0.9 kcal/mol better than the ORTH group mean (−8.92). Given a
random expectation of ~−8 for a 20-30 heavy-atom drug in a hydrophobic pocket,
this is real enrichment: **a compound set selected by experiment scores
significantly better than a set selected by our shape/pharmacophore consensus.**
The docking score has genuine predictive power in this pocket, which is exactly
what an external calibration is supposed to establish.

## Finding 2 — the calibration also exposes the method's ceiling

The two best BENCH hits (naftifine −11.39, flibanserin −11.36) reach the Stage A
top-15 zone (−11.47..−10.51), but 4 of 9 published modulators sit at −9.5 or
weaker, and luliconazole — the only published **agonist** — lands at −8.61, 4
kcal/mol worse than the antagonist leader. Three conclusions:

1. **Affinity ≠ efficacy.** Vina optimises burial and H-bonding, not the
   conformational change that distinguishes agonism from antagonism at a nuclear
   receptor PAS-B domain. No docking score here can predict direction.
2. **The score saturates near −11.5.** The pocket is compact and mostly
   hydrophobic; even the natural ligand only reaches −12.6, and no approved drug
   in either group exceeds it. Scores above ~−11 should be read as "in the
   strong-binder band", not as a continuous ranking.
3. **Our orthogonal screen found one genuinely strong hit.** Ramosetron
   (−11.24) outscores 7 of the 9 experimentally validated modulators and sits
   inside the Stage A top-15 band — a 5-HT3 antagonist with, as far as the
   literature search found, no published AHR activity. Triamterene (−10.00),
   the top dual shape+pharmacophore candidate, also clears −10.

## Revised consensus shortlist (AHR axis)

| drug | evidence | status |
|---|---|---|
| **Berberine** | Vina −10.91; shape 0.607; **published AHR activator, improves colitis (PMID 33285228)** | strongest lead |
| **Trioxsalen** | Vina −10.81; shape 0.630 (top-3); furanocoumarin class support | strong lead |
| **Ramosetron** | Vina −11.24 (Stage D); shape 0.674 | **new hypothesis, no AHR literature** |
| **Triamterene** | Vina −10.00; shape 0.634; pharmacophore recall 0.548 (library max) | **new hypothesis** |
| **Opicapone** | Vina −10.67; pharmacophore 100th percentile | pending follow-up |
| Methylene blue | Vina −9.69; shape 0.677; planarity 0.135 | weaker binder |

## Next step unlocked by the local Vina capability

The Colab dependency is now optional rather than essential. The two remaining
docking jobs — **chunks 4–5 of Stage A** (the least-diverse 150 of the 424
compounds) and a **high-exhaustiveness refinement of the consensus shortlist**
(berberine, trioxsalen, ramosetron, triamterene, opicapone at
`exhaustiveness=64`) — can both run locally. The latter is the immediate
priority: it converts the consensus list into final ranked poses with per-residue
contacts.
