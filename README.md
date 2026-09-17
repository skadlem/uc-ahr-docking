# Project Status — Week 1 of the 90-day plan

*Progress on "start the plan + extend the search." All artifacts below are in
this directory and persist across reboots.*

## Done

### 1. Research campaign extended (12 → 17 angles)
Five new targeted searches, run after the triage pointed at AHR / IL-6:
organoid screening models, **AHR + tryptophan in colitis**, IL-6 blockade,
computational repurposing, and barrier-repair screens. All persisted in
`research/raw_searches/s13–s17`, indexed in `research/INDEX.md`.

### 2. Target triage executed (real result) — `research/TRIAGE_round1.md`
26 genes scored against UC (MONDO_0005101) via the live Open Targets GraphQL API.

**Two findings:**
- **The method revealed its own blind spot.** EGLN1 (ISM012-042 target) scores
  **0.00 genetic** and PDE10A (BEN-8744 target) **0.00 genetic** — both are real
  clinical-stage UC targets found by AI. Pure genetic ranking would have killed
  both. → triage must be **genetics + expression + tractability**, not genetics
  alone. Genetics *did* correctly recover IL23R (0.90) and IL12B (0.88).
- **AHR is the lead**: 0.43 genetic, zero known-drug score, sits on the
  host–microbiome (tryptophan-metabolite) axis, ligand-binding nuclear receptor
  = maximally druggable, not crowded.

### 3. AHR validated by literature — `research/raw_searches/s14_ahr.txt`
- *PNAS* 2020 — microbial tryptophan metabolites regulate **barrier function via AHR**
- ***J Med Chem* 2022** — AHR microbial-metabolite **mimics** alleviate experimental colitis (medchem precedent)
- *Cells* 2022 — AHR therapeutic landscape & safety roadmap
- **IL-6 demoted**: a published case of *exacerbated* UC after anti-IL-6R therapy → riskier lead.

### 4. Structure secured and processed (Weeks 5–6 work, pulled forward)
Queried RCSB programmatically via UniProt → 7 AHR structures. Best docking
starting points:

| PDB | Method | Res | Content |
|---|---|---|---|
| **7ZUB** | cryo-EM | 2.85 Å | AHR + **indirubin co-bound** (chain D) |
| 8QMO | cryo-EM | 2.76 Å | AHR + benzo[a]pyrene co-bound |
| 5NJ8 | X-ray | 3.3 Å | AHR:ARNT active complex |

Extracted to a **docking-ready receptor** (`structures/processed/`):
`AHR_D_receptor.pdb` (2554 atoms, PAS-B domain residues 271–427) +
`indirubin_JY6.pdb` (30 atoms, C₁₆H₁₀N₂O₂ ✓). A co-bound reference ligand means
the docking setup can be **validated by redocking** before any screening.

### 5. Docking pipeline staged
- `scripts/extract_ahr.py` — structure processing ✓ run
- `scripts/redock_validate.sh` — GNINA redocking gate ✓ written
- `scripts/AHR_docking_colab.ipynb` — **one-click GPU notebook** ✓ valid JSON
- `scripts/ot_triage.py`, `scripts/ot_lookup.py` — reproducible triage
- Prep chemistry verified locally (obabel)

### 6. Environment
`.venv` (Python 3.12) with numpy/pandas/scipy/matplotlib/requests/biopython/
rdkit/openbabel — working and verified.

### 7. Notebook de-risked — two real bugs found and fixed

Before handing over the Colab notebook I verified every locally-testable step.
Two bugs in my first draft would have cost you a failed run:

1. **Wrong GNINA download URL.** I had guessed a `v1.3.1/gnina-1.3.1-cuda11.8-…`
   asset name that does not exist. The real assets (checked via the GitHub
   releases API) are `gnina1.3.1` (v1.3.1), `gnina.1.3.2` /
   `gnina.1.3.2.cuda12.8`, and `gnina.cuda12.8.static` (**v1.3.3**). The
   notebook now uses `v1.3.3/gnina.cuda12.8.static`.

2. **Wrong RMSD metric.** My first version used RDKit `GetBestRMS`, which
   *realigns* the two poses — I proved with a self-test that it reports 0.31 Å
   for a pose displaced 1.68 Å from crystal. That would have made bad poses
   look good and defeated the entire validation gate. Replaced with
   `scripts/docking_rmsd.py`, a receptor-frame (non-realigned) heavy-atom RMSD
   with a canonical atom mapping. Self-test: identical = 0.000 Å, translated
   pose = 1.682 Å (exact). Gate: < 2.0 Å = PASS.

### 8. AHR binding pocket characterized — `research/AHR_pocket.md`

Computed all residues within 5 Å of the co-bound indirubin: **25 residues**,
overwhelmingly hydrophobic (Leu×5, Ile×3, Val×3, Phe×4, Ala×2, Pro) with
**His337 closest at 2.8 Å** and **Gln383 / Tyr322 / Ser336** the only plausible
H-bond partners. This is textbook nuclear-receptor LBD architecture —
independent confirmation that AHR is highly druggable — and it gives concrete
pose-inspection criteria for the eventual screen.

### 9. Local GNINA run: hit a real wall (documented, not solved)

I downloaded the CPU-labelled `gnina1.3.1` binary (1.4 GB) to validate the exact
GNINA invocation locally. **It fails to load on this box**: every GNINA release
   binary — even the ones not named "cuda" — links against `libcudart.so.12`,
   `libcusparse.so.12`, `libcublas.so.12` etc. This CPU-only machine has no CUDA
   toolkit, and installing the full CUDA runtime stack via pip repeatedly timed
   out on large wheels. **Colab has CUDA, so this is a local-only limitation,**
   which is exactly why Colab was the right pick.

Consequence: the **only unverified step was the GNINA command itself** — since
verified on Colab (§10: redocking RMSD 0.49 Å). Every upstream step (structure
fetch, chain/ligand extraction, obabel protonation, RMSD math) was already run
and checked locally.

---

## How to run (Colab)

- Validation: `scripts/AHR_docking_colab.ipynb` (gate result already in
  `research/AHR_redocking_result.md`).
- Screen: `scripts/AHR_screening_colab.ipynb` — also at
  https://colab.research.google.com/github/skadlem/uc-ahr-docking/blob/master/AHR_screening_colab.ipynb
  (open → "Copy to Drive" → Runtime → T4 GPU → Run all).

### 10. Redocking GATE PASSED on Colab T4 — `research/AHR_redocking_result.md`

Ran the validation notebook on a Colab T4 (via browser automation, §12):

- **Best-pose RMSD = 0.49 Å** vs crystal indirubin (gate was < 2.0 Å) ✅
- Top CNN score 0.7145; best-pose CNN affinity 6.39
- **Caveat found:** the top-*ranked* pose was a 6.55 Å decoy. CNN scores are
  near-degenerate across binding modes → screening must inspect top 3–5 modes
  per ligand, never just rank #1.

Receptor prep (chain D extraction + obabel protonation), autobox definition, and
GNINA invocation are all trustworthy.

### 11. Screening library built — `screening/library_batch1.sdf`

- 3,417 approved drugs from the **ChEMBL API** (`max_phase=4`) → 2,838 dockable
  (8–55 heavy atoms, no metals, deduplicated)
- Protonated at pH 7.4 (obabel), 3D-embedded with RDKit **ETKDGv3**
- **MaxMin diversity pick of 600 → 589 structures** (2.67 MB)
- Hosted as a **secret GitHub gist** so the notebook stays small:
  `https://gist.githubusercontent.com/skadlem/3b89119ebd49a7cfa0dd5ad6d64ea903/raw/library_batch1.sdf`
- Script: `scripts/build_library.py`

### 12. Colab execution via browser automation (chrome-devtools MCP)

Drove the user's own Chrome (CDP port 9222) through the chrome-devtools MCP with
`--workspace=/home/madiyar/uc`. Working flow (cells run headlessly, results read
back via DOM):

1. Notebook pushed to public repo **github.com/skadlem/uc-ahr-docking**, then
   opened via its Colab–GitHub URL (the file-picker upload path was unreliable;
   URL navigation + "Copy to Drive" is robust).
2. Free-tier Colab allows **one** runtime at a time — the sessions dialog
   ("Завершить другие сеансы") must terminate the previous notebook's runtime
   before a new notebook can connect a T4.
3. Colab's toolbar buttons (Connect T4, Run all) live in **shadow DOM**; the MCP
   a11y click can't reach them — drive them via `evaluate_script` + dispatched
   `MouseEvent` on `colab-connect-button.shadowRoot`, and cells via
   `colab.global.notebookModel.cells[i].setText()` + `cells[i].runButton`.
4. Code cells can be edited programmatically (above), so the protocol iterated
   **without** re-uploading or losing the warm runtime.

### 13. Throughput reality on free Colab T4 (measured, not assumed)

The T4 instance has **only 2 vCPUs**, and GNINA's Monte-Carlo sampling is
CPU-bound — speed is set by ligand size, not the GPU:

| Setting | Small ligands | Mid-library (≤46 HA) |
|---|---|---|
| exh 8 + CNN rescore + 3 modes | — | **40 s/ligand** (→ 6.5 h for 589) |
| exh 4, no CNN, 1 thread | 9.6 s | — |
| exh 2, no CNN, 1 thread | 5.3 s | **35 s/ligand** |

Consequences that shaped the screen:

- `--cpu 1` + `OMP_NUM_THREADS=1` is **faster** than `--cpu 4` on 2 vCPUs
  (thread contention); 2 such processes in parallel saturate the box.
- The naive "screen all 589 with CNN" plan would have taken 6+ h with **no
  checkpointing** (GNINA writes its output SDF only at the end) — unacceptable.
- → the screen is **hierarchical** (§14).

### 14. Hierarchical screen — running

`scripts/make_hierarchical_notebook.py` → `scripts/AHR_screening_colab.ipynb`:

- **Stage A** — Vina-only breadth pass (`--cnn_scoring none --exhaustiveness 2`,
  `--autobox_add 5`, `--cpu 1`), library filtered to **≤ 34 heavy atoms**
  (**424/589 pass** — the PAS-B pocket is compact, indirubin is 20 HA, so larger
  molecules rarely fit and dominate runtime), chunked 6 × ~75 with **2 parallel
  workers → per-chunk checkpoints** (`stageA_*.sdf.gz`) and live ETA.
- **Stage B** — CNN refinement (`--cnn fast --cnn_scoring rescore`,
  `--exhaustiveness 12`, `--num_modes 3`) of the **top 30** Vina hits →
  `refined.sdf.gz`.
- **Stage C** — rank by CNN affinity + pocket contact vs His337 / Gln383 /
  Tyr322 → `screen_results.csv`.

## Suggested next steps after the screen finishes
- Triage Stage A hits: check the top hits for AHR literature support, then
  re-dock the best 10 at `--exhaustiveness 64` as a final pose check.
- Stage 2 triage: add UC single-cell expression (GSE125527) + tractability to the
  composite score; re-rank.
- Mendelian randomization on the top 5 targets for causal direction.
- REINVENT 4 de novo run with docking score in the scoring function.
