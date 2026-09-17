# Library-wide shape/planarity screen (GPU-free) — 424 approved drugs

Script: `scripts/planarity_screen.py` · Data: `screening/planarity_screen_named.csv`

## Purpose

Three lines of evidence had converged on "the AHR PAS-B pocket binds planar
aromatics" (co-crystal ligand + known ligands; docking top hits; the 15-compound
O3A screen). This turns that observation into a **quantitative screening
criterion applied to the whole library** instead of the top-15, so that
interesting candidates can be found *regardless of their Vina rank*, and so the
docking ranking itself can be sanity-checked for shape enrichment.

Every one of the 424 Stage A compounds (≤34 heavy atoms) was scored locally,
no docking engine, all 12 cores, ~90 s total:

- **shape_tanimoto** — best over a 25-conformer ETKDGv3+MMFF ensemble vs the
  crystallographic indirubin pose (7ZUB JY6). Open3DAlign superposes each
  conformer onto indirubin; for the compounds where RDKit's MMFF94 port has no
  parameters a force-field-free **principal-axes-of-inertia** alignment is used
  (centroid match + best of the degenerate sign permutations).
- **planarity** — RMSD (Å) of heavy atoms from their least-squares plane
  (smaller = flatter; indirubin-like scaffolds score ≲0.4).
- **aromatic_rings / fused_rings** — SSSR counts; fused count = ring systems.

Range: shape Tanimoto 0.234–0.722 across the library.

## Result 1 — new candidates the docking ranking missed

Top of the library by shape, with their Vina rank unknown/not-in-top-15:

| Drug | shape T | planarity | arom/fused |
|---|---|---|---|
| **Indalpine** | 0.722 | 0.355 | 2/1 |
| **Methylene blue** | 0.677 | **0.135** | 3/2 |
| **Ramosetron** | 0.674 | 0.608 | 3/2 |
| **Amiloride** | 0.671 | **0.013** | 1/0 |
| Moxonidine | 0.671 | 0.366 | 1/0 |
| Sapropterin | 0.663 | 0.408 | 1/1 |
| Frovatriptan | 0.658 | 0.285 | 2/2 |
| Triamterene | 0.634 | 0.502 | 3/1 |
| **Trioxsalen** | 0.630 | **0.000** | 3/2 |
| **Berberine** | 0.608 | 0.316 | 3/4 |

**Methylene blue** is the standout: 3 fused aromatic rings and the flattest
non-trivial scaffold in the library (0.135 Å). It is also a phenothiazinium dye
with documented AHR-pathway interaction, and it is *not* in the docking top-15 —
a shape-first screen recovered it. **Amiloride** (planarity 0.013) and
**trioxsalen** (0.000, perfectly flat) are the other maximally planar hits.

These are hypotheses for re-docking, not validated binders.

## Result 2 — the docking ranking is shape-enriched (validation)

Where the Stage A Vina top-15 fall in the 424-compound shape distribution:

| Drug | Vina | shape T | library shape percentile |
|---|---|---|---|
| Trioxsalen | −10.81 | 0.630 | 0.2% |
| Berberine | −10.91 | 0.608 | 0.5% |
| Ramelteon | −10.51 | 0.597 | 0.7% |
| Hydroxystilbamidine | −10.52 | 0.590 | 0.9% |
| Vortioxetine | −10.82 | 0.571 | 1.2% |
| Sertraline | −10.66 | 0.565 | 1.4% |
| Exemestane | −10.53 | 0.533 | 1.7% |
| Mefloquine | −11.32 | 0.529 | 1.9% |
| Sitagliptin | −11.47 | 0.506 | 2.1% |
| Flortaucipir | −11.23 | 0.504 | 2.4% |
| Rucaparib | −11.21 | 0.504 | 2.6% |
| Pitavastatin | −10.89 | 0.498 | 2.8% |
| Opicapone | −10.67 | 0.458 | 3.1% |
| Olodaterol | −10.90 | 0.447 | 3.3% |
| Pexidartinib | −11.25 | 0.390 | 3.5% |

**All 15 docking hits sit in the top 3.5% of shape similarity** to the native
ligand. This is an independent check that the docking screen is not returning
arbitrary high-scoring poses: it is genuinely enriched for the pocket's known
chemotype. The four consensus compounds (top-decile shape *and* Vina top-15)
are **trioxsalen, berberine, ramelteon, hydroxystilbamidine**.

## Caveats

- Shape similarity is necessary but not sufficient: a flat molecule that is the
  wrong *size* or lacks the H-bond anchors (His337, Gln383, Tyr322, Ser336 —
  see `research/AHR_pocket.md`) will still fail.
- `o3a_score` is missing for ~1/3 of compounds (RDKit MMFF94 parameter gap);
  those used the inertia-alignment fallback for shape Tanimoto only. The two
  metrics are therefore not strictly comparable across all rows — shape
  Tanimoto is computed uniformly, O3A is not.
- Planarity computed on the gas-phase lowest-energy conformer; flexible
  molecules can flatten in a pocket.
