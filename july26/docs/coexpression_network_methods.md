# Co-expression network & guilt-by-association candidate discovery

## Objective
Build a co-expression network across the entire RNA-seq dataset, then use the
network to nominate candidate genes **outside** the curated genes-of-interest
(GOI) list that may participate in the same light-response biology, ranked by
how likely they are involved in each treatment.

## Input
- Expression: `data/matrices/vst_matrix.csv` — DESeq2 variance-stabilized values,
  12,723 expressed genes × 12 libraries (all four conditions, 3 reps each).
- GOI: 257 publication-derived FGSG_ genes; **245** present in the network
  (12 dropped at the DESeq2 expression filter due to near-zero counts:
  FGSG_00056, FGSG_01776, FGSG_02395, FGSG_03334, FGSG_03381, FGSG_06467,
  FGSG_09090, FGSG_09132, FGSG_10458, FGSG_10460, FGSG_10463, FGSG_12127).
- Treatment DEG sets (vs control, DESeq2, padj<0.05 & |log2FC|≥1), oriented so
  positive = up in treatment: black (39 GOI DEG), dark (83), low (85).

## Network construction (WGCNA methodology, Python/NumPy)
`bioconductor-wgcna`/`r-wgcna` were not installable in the project R 4.5.3
environment, so the standard WGCNA workflow was reproduced directly:
1. **Correlation** — Pearson correlation across the 12 samples for all pairs.
2. **Signed adjacency** — `a_ij = ((1 + r_ij) / 2)^β`.
3. **Soft-threshold selection** — scanned β = 1..30; signed scale-free-topology
   R² crosses 0.8 at β=14 and plateaus ~0.85–0.88. Chose **β = 18**
   (R² = 0.849, mean connectivity ≈ 201.6): above threshold and WGCNA's
   recommended power for signed networks at n < 20.
4. **TOM** — Topological Overlap Matrix from the adjacency (off-diagonal range
   ~2.3e-10 to 0.8135).
5. **Modules** — average-linkage hierarchical clustering on `1 − TOM`,
   `fcluster` at height 0.985, clusters < 30 genes merged into grey/unassigned.
   Result: **15 modules + 450 unassigned**. Modules 1 (4,891 genes) and 2
   (4,714) dominate, consistent with PC1 capturing 80% of variance.

## Candidate extraction
GOI×non-GOI TOM block = 245 × 12,478. Applied a stringent overlap cut of
**TOM ≥ 0.40** (~99th percentile). Result: **32,730 GOI-anchored edges** →
**2,982 candidate genes** (each with ≥1 strong edge to a GOI).

## Per-treatment ranking (guilt-by-association)
For each treatment t and candidate c:
- `gba_weight(c,t) = Σ TOM(c, g)` over GOI g that are DEG in treatment t.
- `score(c,t) = 0.6·rank_pct(gba_weight) + 0.4·rank_pct(|own_log2FC| · own_sig)`
  where own_sig is 1 if c itself is a DEG in t (padj<0.05 & |log2FC|≥1), else 0.
- Combined score = mean of the three per-treatment scores;
  `n_treatments_sig` counts treatments where c is itself a DEG.

**1,070 candidates are DE in ≥1 treatment; 115 in all three.**

### Top candidate per treatment (by GBA weight)
| Treatment | Top candidate | GBA weight | own log2FC |
|-----------|---------------|-----------:|-----------:|
| Black light vs control | FGSG_01379 | 6.23 | −4.69 |
| Dark light vs control  | FGSG_07721 | 16.11 | −6.44 |
| Low light vs control   | FGSG_09697 | 20.16 | −6.26 |

### Highest combined-ranked candidates (all DE in all 3 treatments)
FGSG_08049, FGSG_07530 (alcohol oxidase), FGSG_01823, FGSG_11888, FGSG_06468,
FGSG_11270, FGSG_03306, FGSG_07921, FGSG_05348, FGSG_03638, FGSG_04971,
FGSG_01379. Most reside in module 2 and connect to a shared cluster of GOI
(the FGSG_023xx and FGSG_036xx neighborhoods).

## Interpretation notes
- The two dominant modules reflect the strong global light-response axis (PC1).
  TOM values are correspondingly inflated, which is why a high absolute cut
  (0.40) rather than a fixed nominal WGCNA threshold was used.
- Candidates are nominations for follow-up, not validated members of any
  pathway. The GBA logic assumes co-expression neighbors of treatment-responsive
  GOI are themselves plausible participants in that treatment's response.

## Outputs
- Tables: `results/coexpression_network/` (gene_modules, edges, candidate_genes,
  candidate_ranking_{combined,black,dark,low}, goi_treatment_membership).
- Figures: `figures/network_wgcna_diagnostics.png`,
  `figures/candidate_ranking_barplots.png`,
  `figures/coexpression_network_goi_candidates.png`.

## Candidate functional annotation (GO over-representation)

To interpret the top candidates (most annotated only as "uncharacterized
protein"), a hypergeometric GO over-representation test was run on the **top-40
combined-ranked candidates** (33 GO-annotated), against a universe of the 8,961
GO-annotated genes, with Benjamini-Hochberg FDR correction. Terms with ≥2
candidate hits were tested (13 terms); **6 are significant at padj < 0.05**:
phosphatidylcholine lysophospholipase A1 activity / phospholipid catabolic
process (FGSG_02823, FGSG_02824), iron-sulfur cluster binding (FGSG_03026,
FGSG_08425), FAD binding (FGSG_03638, FGSG_07530, FGSG_13840), transmembrane
transport (7 genes), and integral component of membrane (17 genes). Results:
`results/coexpression_network/candidate_top40_GO_enrichment.csv` (each row lists
the candidate FGSG IDs carrying that term).

**Per-term subnetworks.** For each tested GO term an induced subgraph was built
from the annotated candidates plus all their TOM ≥ 0.40 neighbors, rendered as a
13-panel figure (`figures/candidate_GO_subnetworks_13.png`) with GOI colored by
functional category and the term-carrying candidates ringed in gold. Panels that
separate into two components indicate the candidates span both dominant WGCNA
modules (independent co-regulation programs). Full membership:
`results/coexpression_network/candidate_GO_subnetwork_membership.csv`.

**GOI functional categories.** `results/coexpression_network/goi_category_map.csv`
assigns each of the 257 GOI to one of six categories (Pathogenicity,
Photoreceptors, Pigment, Toxins, Metabolism, Effectors) using the curated v2
list priority, reproducing `GOI_overlap_by_category.csv`.
