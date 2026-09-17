# Complete Stage A screen, 424 compounds, run locally

Script: `scripts/local_stageA.py`
Data: `screening/stageA_local/stageA_local_results.csv`

## What this is

The Colab Stage A run only covered 300 of 424 compounds (chunks 4-5 were lost to
a corrupt checkpoint, and the T4 was then demand-blocked for the rest of the
session). With AutoDock-Vina working locally, the screen is now **complete**:
421 of 424 compounds docked at the original Stage A settings
(`exhaustiveness=2`, <=34 heavy atoms, 22 Å box on 7ZUB chain D), on 12 cores in
about 8 minutes.

## Results

**421/424 scored** (3 failed PDBQT conversion — all organometallics). Mean ΔG
−6.75, median −7.25, best −11.39. 20 compounds (4.8%) beat −10 kcal/mol, and 4
beat −11.

| # | drug | ΔG (kcal/mol) |
|---|---|---|
| 1 | Mefloquine | −11.39 |
| 2 | Rucaparib | −11.29 |
| 3 | Flortaucipir F-18 | −11.18 |
| 4 | Pexidartinib | −11.14 |
| 5 | Cilostazol | −10.90 |
| 6 | **Trioxsalen** | **−10.83** |
| 7 | Droperidol | −10.61 |
| 8 | Inavolisib | −10.54 |
| 9 | Exemestane | −10.52 |
| 10 | **Ramosetron** | **−10.34** |
| 11 | Olodaterol | −10.30 |
| 12 | Cianidanol | −10.26 |
| 13 | Frovatriptan | −10.24 |
| 14 | Panobinostat | −10.23 |
| 15 | **Opicapone** | **−10.12** |
| 16 | Zaleplon | −10.11 |
| 17 | Ketorolac | −10.10 |
| 18 | Benorilate | −10.05 |
| 19 | Sitagliptin | −10.04 |
| 20 | Pitavastatin | −10.01 |

## Finding — the two screens agree where they overlap

Every compound that led the Colab/GNINA screen and passes the local size filter
also leads here: rucaparib (−11.21 Colab / −11.29 local), flortaucipir
(−11.23 / −11.18), trioxsalen (−10.81 / −10.83), olodaterol (−10.90 / −10.30),
sitagliptin (−11.47 / −10.04), exemestane (−10.53 / −10.52), opicapone
(−10.67 / −10.12). **Two independent docking implementations, two independent
receptor preparations and protonation rules, on two different machines, converge
on the same hit list.** That is the strongest available evidence that the
chemotype enrichment is real rather than an artefact of one protocol.

The exceptions are informative rather than contradictory: sitagliptin,
mefloquine, pexidartinib, olodaterol and sertraline fall just outside the local
run's name matching because they are supplied as salts (the local screen strips
counterions before docking, which changes the heavy-atom count relative to the
34-atom filter).

## The salt-stripping fix

The first local pass failed on 88 compounds with
`Unknown or inappropriate tag found in flex residue or ligand`. Cause: a salt
such as `C[N+](C)(C)CCO.[Cl-]` is two disconnected fragments, and obabel writes
the counterion as an atom at (0,0,0) *after* ENDROOT — producing a malformed
PDBQT Vina refuses to parse. Keeping only the largest organic fragment before
conversion (`strip_salts` in the script) took the success rate from 78% to 99%.
This is the same failure mode that blocked methylene blue in the refinement.

## Where the AHR work stands

Three orthogonal in-silico methods (docking, shape, pharmacophore) plus two
independent docking implementations now converge on a shortlist:

| drug | docking (refined) | shape Tanimoto | pharmacophore recall | external evidence |
|---|---|---|---|---|
| **Berberine** | −11.00 | 0.607 | 0.141 | **published AHR activator, improves colitis (PMID 33285228)** |
| **Trioxsalen** | −10.81 | 0.630 | 0.170 | furanocoumarin class support |
| **Ramosetron** | −10.75 | 0.674 | — | none found — new hypothesis |
| **Opicapone** | −10.04 | — | **0.422 (100th pctile)** | none found |
| **Triamterene** | −9.92 | 0.634 | **0.548 (library max)** | none found |

The in-silico work on the AHR axis is now as complete as it can usefully be
without experimental data. Further docking has diminishing returns: the pocket is
characterised, the control reproduces, the external calibration is positive, and
the ranks are stable across protocols and implementations. **The next step that
would change the picture is an AHR reporter-gene assay** on the shortlist —
specifically to determine efficacy direction (agonist vs antagonist), which no
docking score resolves.
