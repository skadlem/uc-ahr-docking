# Round-2 target triage: genetics + small-molecule tractability

**Script:** `scripts/ot_triage_stage2.py` -> `research/TRIAGE_round2.csv`
**Source:** Open Targets GraphQL v4, UC = `MONDO_0005101`

Round 1 ranked targets on **genetic association alone** and exposed the method's
blind spot: EGLN1 and PDE10A - both real clinical-stage, AI-discovered UC targets
- score **0.000** genetically and would have been discarded. Pure genetic ranking
kills druggable, non-genetic mechanisms. This round adds the axis that rescues
them: **small-molecule tractability**.

## Composite = 0.45 x genetic + 0.20 x overall association + 0.35 x tractability(SM)

| rank | target | genetic | overall | tract(SM) | composite | SM evidence |
|------|--------|---------|---------|-----------|-----------|-------------|
| 1 | JAK2 | 0.818 | 0.723 | 1.0 | **0.863** | Approved Drug |
| 2 | IL12B | 0.883 | 0.732 | 0.6 | **0.754** | Structure with Ligand |
| 3 | IL10 | 0.883 | 0.548 | 0.4 | **0.647** | Med-Quality Pocket |
| 4 | STAT3 | 0.702 | 0.437 | 0.6 | **0.613** | HQ Ligand + Structure w/ Ligand |
| 5 | **AHR** | 0.432 | 0.267 | **1.0** | **0.598** | **Approved Drug** |
| 6 | IL6R | 0.578 | 0.356 | 0.6 | 0.541 | Structure with Ligand |
| 7 | IL23R | 0.904 | 0.558 | **0.0** | 0.518 | none (biologic-only) |
| 8 | TYK2 | 0.052 | 0.588 | 1.0 | 0.491 | Approved Drug |
| 9 | PDE10A | 0.000 | 0.104 | 1.0 | 0.371 | Approved Drug |
| 10 | EGLN1 | 0.000 | 0.004 | 1.0 | 0.351 | Approved Drug |
| 11 | GPX4 | 0.035 | 0.025 | 0.6 | 0.231 | HQ Ligand |
| 12 | HIF1A | 0.000 | 0.042 | 0.6 | 0.218 | HQ Ligand |

## What this changes

**AHR is the best *unexploited small-molecule* target in the set.** Every target
that outranks it is either already a marketed UC drug target (JAK2 = tofacitinib
/ filgotinib; IL12B+IL23R = ustekinumab), a cytokine that cannot be agonised with
a small molecule (IL10), or a transcription factor with no direct-binding pocket
(STAT3). AHR is the only top-5 target that combines:

- **maximal small-molecule tractability** (Open Targets "Approved Drug" SM tier,
  plus High-Quality Ligand, Structure with Ligand, Druggable Family),
- **no approved UC drug** acting on it,
- an **experimentally validated binding pocket** - our GNINA redocking to
  crystallographic indirubin in 7ZUB reproduced the pose to **0.49 A**.

**IL23R illustrates the converse**: the strongest genetics in the set (0.904) but
zero small-molecule tractability - it is a cell-surface receptor reachable only
by a biologic. Genetics alone would have ranked it #1 for a small-molecule
campaign, which is exactly the trap this round was designed to catch.

**The zero-genetic controls survive**: EGLN1 and PDE10A remain at 0.000 genetic
yet sit at the top SM tractability tier with approved drugs. The tractability
axis is what keeps them visible.

## Methodological notes / bugs caught

1. **`target.tissueExpression` was removed** from the API; its replacement
   `baselineExpression{rows{tissueBiosample{biosampleName} median}}` returns
   per-cell-type rows (AHR: 25 rows, mostly pancreas/pituitary) whose tissue
   hierarchy does not roll up to gut in one call. The planned gut-expression axis
   was **dropped rather than reported unreliably**; every candidate is gut-
   expressed by construction of the round-1 shortlist, so the axis was
   low-variance.
2. **Tractability schema changed**: `tractability{ smallMolecule{ value } }` ->
   a flat `tractability{ label modality value }` list of Boolean buckets
   (`modality` in SM/AB/PR/OC), collapsed here into the tiered score above.
3. **A hardcoded Ensembl id silently returned the wrong gene**: `ENSG00000156234`
   is **CXCL13**, not IL12B - it scored 0.024 and would have buried a top target.
   Correct IL12B = `ENSG00000113302` (IL10 = `ENSG00000136634`). IDs are now
   **resolved by symbol at runtime** and cross-checked against seeds, so a stale
   id fails loudly instead of corrupting the ranking.
