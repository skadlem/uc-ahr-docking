# Pharmacophore screen (Gobbi 2D) — feature-level validation of the chemotype

Script: `scripts/pharmacophore_screen.py` · Data: `screening/pharmacophore_screen.csv`

## Why

Shape asks "does it fill the same space?"; pharmacophore asks "does it present
the same *functional features* — donors, acceptors, hydrophobic centres, aromatic
rings, cations, anions — in the same arrangement?" A flat greasy molecule can
match indirubin's shape while carrying none of its hydrogen-bond machinery, so a
credible prioritisation needs both axes.

Method: **Gobbi 2D pharmacophore fingerprints** (Gobbi & Lee, *JCIM* 2003,
doi:10.1021/ci0256062) as implemented in RDKit `Pharm2D`. Each molecule is encoded
as a bitstring over all feature triples with pairwise topological distances
binned (0-2, 2-3, ... 8+ bonds). Similarity to indirubin = Tanimoto over the
39,972-bit space. Because distances are topological, this is
conformation-independent and therefore *independent of* the 3D shape metric —
genuinely orthogonal information.

Two numbers are reported:
- **pharm_tanimoto** — full-space Tanimoto (intrinsically low across such a large
  bit space; median 0.051, max 0.221 — useful for ranking, not as an absolute)
- **pharm_recall** — of indirubin's 135 on-bits, what fraction does this compound
  present *at all*. This is the interpretable one: it asks whether the H-bond
  donor/acceptor pattern the pocket actually uses is present.

## Result 1 — the metric separates hydrophobic-burial false positives

Stage A Vina hits, ranked by pharmacophore recall percentile within the library:

| Drug | Vina ΔG | pharm_recall | percentile |
|---|---|---|---|
| Opicapone | −10.67 | **0.422** | 100% |
| Rucaparib | −11.21 | **0.393** | 93% |
| Flortaucipir F-18 | −11.23 | 0.363 | 87% |
| Pexidartinib | −11.25 | 0.296 | 80% |
| Olodaterol | −10.90 | 0.289 | 73% |
| Sitagliptin | −11.47 | 0.267 | 63% |
| Mefloquine | −11.32 | 0.267 | 63% |
| Pitavastatin | −10.89 | 0.200 | 53% |
| Trioxsalen | −10.81 | 0.170 | 47% |
| Berberine | −10.91 | 0.141 | 40% |
| Hydroxystilbamidine | −10.52 | 0.126 | 33% |
| Vortioxetine | −10.82 | 0.104 | 23% |
| Ramelteon | −10.51 | 0.104 | 23% |
| Sertraline | −10.66 | 0.052 | 13% |
| Exemestane | −10.53 | **0.000** | 7% |

**Interpretation.** The docking score rewards hydrophobic burial, which is
correct (the pocket is mostly Leu/Ile/Val/Phe) but incomplete. Exemestane, the
14th-ranked docking hit, presents **zero** of indirubin's pharmacophore features —
it is a greasy steroid filling volume. Sertraline, ramelteon and vortioxetine are
the same story. Meanwhile opicapone and rucaparib — which carry lactam/nitrile
donor-acceptor arrays matching indirubin's own — sit at the 93rd-100th percentile.
**Pharmacophore recall is the missing filter on a hydrophobic-pocket docking
screen**, and it cheaply demotes the docking hits that are least likely to be
specific binders.

Note the honest counter-example: **berberine** has only 40th-percentile recall
(yet is the experimentally validated AHR activator). It is a quaternary
protoberberine — no N-H, only methoxy/methylenedioxy oxygens — so it *cannot*
reproduce indirubin's lactam donor pattern and activates AHR without it. The
pocket is hydrophobic-dominated; pharmacophore recall is a filter, not a
necessity.

## Result 2 — consensus ranking with shape

`consensus = 0.45·shape_rank + 0.30·pharm_tanimoto_rank + 0.25·pharm_recall_rank`
(shape weighted highest because it carries the pocket-volume evidence).

Top of the 421-compound consensus:

| Drug | consensus | shape T | pharm_recall | planarity |
|---|---|---|---|---|
| **Triamterene** | 0.949 | 0.634 | **0.548** | 0.502 |
| **Methylene blue** | 0.912 | 0.677 | 0.200 | 0.135 |
| **Amiloride** | 0.901 | 0.671 | 0.348 | 0.013 |
| Tipiracil | 0.897 | 0.592 | 0.267 | 0.542 |
| Acyclovir | 0.895 | 0.622 | 0.356 | 0.359 |
| Minoxidil | 0.891 | 0.646 | 0.319 | 0.397 |
| Methaqualone | 0.891 | 0.605 | 0.207 | 0.743 |
| Isoxicam | 0.891 | 0.612 | 0.356 | 0.502 |
| Lenalidomide | 0.889 | 0.604 | 0.296 | 0.720 |
| Olanzapine | 0.861 | 0.564 | 0.393 | 0.657 |

**Triamterene** is the standout: 3 aromatic rings, 0.634 shape similarity, and
the **highest pharmacophore recall in the entire library (0.548)** — it presents
more of indirubin's functional pattern than any other approved drug screened. It
was nowhere near the docking top-15. These are hypotheses for re-docking
(Stage D), not validated binders.

## Limitations

- 2D topological distances are a proxy for 3D inter-feature geometry; a molecule
  can present the right feature *set* in the wrong spatial arrangement.
- ~1/3 of compounds have no RDKit MMFF94 parameters (so no O3A colour score);
  shape Tanimoto is uniform across all compounds but O3A is not.
- Consensus weights are a judgement call, not a fitted parameter; with no
  confirmed-activity labels for the orthogonal hits they cannot be optimised yet.
