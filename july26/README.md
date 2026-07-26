# *Fusarium graminearum* PH-1 — Light-Regime RNA-seq Differential Expression

RNA-seq processing and differential expression for *Fusarium graminearum* PH-1
grown under four light regimes, mapped to the **NCBI RefSeq** reference genome
and annotation (assembly `GCF_000240135.3`, ASM24013v3) using native
**FGSG_ gene IDs**. Differential expression is reported three ways — from raw
read counts (DESeq2) and from TPM and FPKM (limma-trend) — for **all six
pairwise contrasts**.

## Experimental design

| Condition | Code | Description | Replicates |
|-----------|------|-------------|------------|
| Control   | H    | Half/half — 12 h light / 12 h dark (**reference level**) | H13_S7, H14_S8, H15_S9 |
| Black light | B  | Black light | B10_S1, B11_S2, B12_S3 |
| Dark light | D   | Dark light | D10_S4, D11_S5, D12_S6 |
| Low light | L    | Low light | L13_S10, L14_S11, L15_S12 |

12 paired-end libraries (~101 bp), 3 biological replicates per condition.

## Pipeline

```
FASTQ
  │  fastp 1.1.0            adapter/quality trim (Q20, len≥36, --detect_adapter_for_pe)
  ▼
trimmed reads
  │  HISAT2 2.2.2           --rna-strandness RF --dta, spliced align to NCBI genome
  ▼
sorted BAMs (samtools 1.22.1)
  ├─ featureCounts 2.1.1    -p --countReadPairs -s 2 -t exon -g gene_id  → raw counts
  └─ StringTie 3.0.3        -e -B --rf                                    → TPM, FPKM
  ▼
gene × sample matrices (13,725 gene models × 12 libraries, FGSG_ IDs)
  ├─ DESeq2 (R 4.5.3)       raw counts, Wald test + apeglm shrinkage
  ├─ limma-trend            log2(TPM+1),  eBayes(trend=TRUE, robust=TRUE)
  └─ limma-trend            log2(FPKM+1), eBayes(trend=TRUE, robust=TRUE)
  ▼
DE tables · VST clustered heatmaps · hypergeometric GO enrichment (blast2go GAF)
```

**Significance threshold (all methods):** adjusted *p* < 0.05 **and** |log2FC| ≥ 1.
Positive log2FC = up in the first-named condition of each contrast.

## Quality control

QC is strong: **PC1 captures 80%** of variance (PC2 10.4%), within-condition
replicate correlation r ≈ 0.99 (0.928–0.993 overall), and samples cluster
cleanly by light regime. Per-library alignment rate ~97.5%. One library,
**B12_S3**, shows 30% multimapping (a known property of that library) and the
smallest DESeq2 size factor (0.648); it was retained, and multimapping reads
were excluded from counting.

![Alignment rates](figures/qc_alignment_rates.png)
![PCA](figures/qc_pca.png)
![Sample correlation](figures/qc_correlation.png)
![Library size and gene detection](figures/qc_libsize_detection.png)

## Differential expression

DEG counts per contrast (padj < 0.05 & |log2FC| ≥ 1):

| Contrast (B vs A) | DESeq2 (counts) | limma (TPM) | limma (FPKM) |
|-------------------|:---------------:|:-----------:|:------------:|
| dark_light vs black_light  | 1,989 | 1,965 | 1,743 |
| control vs black_light     |   753 |   780 |   639 |
| low_light vs black_light   | 2,721 | 2,877 | 2,681 |
| control vs dark_light      | 1,219 | 1,116 |   950 |
| low_light vs dark_light    |   631 |   617 |   526 |
| low_light vs control       | 1,639 | 1,634 | 1,436 |

Across DESeq2 contrasts: **3,525 unique DEGs**; 479 shared across the three
light-vs-black contrasts; **42 DE in all six**. The three quantification
methods agree well — pooled DESeq2-vs-limma(TPM) log2FC Pearson r = **0.94**,
with 48–65% three-way overlap of significant genes per contrast.

![Volcano plots](figures/volcano_counts_grid.png)
![MA plots](figures/ma_counts_grid.png)
![Method concordance](figures/method_concordance.png)

## Expression heatmaps

Clustered heatmaps built from VST-transformed raw counts across **all 12,723
expressed gene models**, row z-scored, Ward hierarchical clustering on both axes.

![All gene models](figures/heatmap_all_genes_clustered.png)
![DEG union](figures/heatmap_DEGs_clustered.png)
![Top-50 variable DEGs](figures/heatmap_top50_DEGs_labeled.png)

## GO enrichment

Hypergeometric over-representation against the OmicsBox blast2go GAF
(20230802), BP/MF/CC tested separately, BH-FDR, background = genes tested in
each contrast. **168 significant GO terms** (padj < 0.05) across the six
contrasts. The dominant enriched theme is oxidation-reduction (P450 /
monooxygenase, heme/iron binding, catalase / hydrogen-peroxide catabolism —
oxidative-stress detoxification) plus amino-acid transport.

![GO dotplot](figures/go_dotplot.png)

## Genes of interest & treatment-vs-control overlap

257 publication-derived genes of interest — all 257 present in the NCBI
reference — were cross-referenced against the DE results; **141 are DEG in ≥1
contrast**. Treatment-vs-control DEG overlap (DESeq2): union 2,369 genes,
**197 DEG in all three treatments** vs control.

![Venn — treatment vs control](figures/venn_treatment_vs_control.png)

## Co-expression network & guilt-by-association candidates

A signed co-expression network was built across the **entire dataset** (all
12,723 expressed genes × 12 libraries) using WGCNA methodology implemented in
Python/NumPy. Pearson correlations were transformed to a signed adjacency
`a = ((1 + r) / 2)^β`; the soft-threshold power **β = 18** was selected from a
scale-free-topology scan (signed R² = 0.85, mean connectivity ≈ 202) — the
recommended power for signed networks at this sample size. A Topological
Overlap Matrix (TOM) was computed and average-linkage clustering on `1 − TOM`
yielded **15 modules** (+450 unassigned); two large modules (4,891 and 4,714
genes) dominate, reflecting the 80% of variance captured by PC1.

![WGCNA network diagnostics](figures/network_wgcna_diagnostics.png)

**Candidate discovery.** GOI-anchored edges were extracted from the 245-of-257
GOI present in the network (12 dropped at the expression filter). Using a
stringent cut of **TOM ≥ 0.40** (~99th percentile of GOI→non-GOI overlap),
**32,730 edges** connect GOI to **2,982 candidate genes** outside the GOI list.

**Per-treatment ranking (guilt-by-association).** For each treatment a
candidate's GBA weight is the summed TOM to the GOI that are DEG in that
treatment; the combined score blends GBA rank (0.6) with the candidate's own
DE magnitude and significance (0.4). **1,070 candidates are DE in ≥1
treatment; 115 in all three.** Top candidate per treatment: black light
**FGSG_01379** (own log2FC −4.7), dark light **FGSG_07721** (−6.4), low light
**FGSG_09697** (−6.3). Highest combined-ranked candidates include
**FGSG_08049** and **FGSG_07530** (alcohol oxidase) — both DE in all three
treatments and hubs of the shared module-2 co-expression neighborhood.

![Candidate ranking per treatment](figures/candidate_ranking_barplots.png)

![GOI-anchored co-expression subnetwork](figures/coexpression_network_goi_candidates.png)

**GOI functional categories.** The subnetwork re-colored by the curated GOI
functional categories (Pathogenicity, Photoreceptors, Pigment, Toxins,
Metabolism, Effectors); candidates are shown as neutral diamonds. Category
assignment reproduces `results/genes_of_interest/GOI_overlap_by_category.csv`.

![GOI subnetwork by functional category](figures/coexpression_network_by_category.png)

**What are the candidates doing? (GO over-representation).** Because most top
candidates are annotated only as "uncharacterized protein," we ran a
hypergeometric GO over-representation test (BH-FDR) on the **top-40 combined
candidates** (33 of 40 GO-annotated; universe = 8,961 GO-annotated genes).
**Six GO terms are significantly enriched (padj < 0.05):** phospholipid
catabolism / phosphatidylcholine lysophospholipase A1 (FGSG_02823/02824),
iron-sulfur cluster binding, FAD binding, transmembrane transport, and integral
component of membrane. The signal points at **redox machinery, membrane
transport, and lipid catabolism** — echoing the dominant oxidation-reduction /
transmembrane-transport signature of the light-response DE. With only 40 genes
this enrichment is hypothesis-generating rather than definitive.

![Category-colored subnetwork with candidate GO enrichment](figures/coexpression_network_by_category_with_GO.png)

**Per-term candidate subnetworks.** For each of the 13 tested GO terms, the
candidates carrying that term (gold-ringed diamonds) are shown together with
every gene they connect to (TOM ≥ 0.40). Panels that split into two islands
(e.g. transmembrane transport, integral component of membrane) do so because
those candidates belong to **both dominant WGCNA modules**, which have no strong
cross-module edges — i.e. two independent co-regulation programs both enriched
for the same function.

![Per-GO-term candidate subnetworks](figures/candidate_GO_subnetworks_13.png)

## Repository layout

### Documentation
- [docs/README_ncbi_analysis.md](docs/README_ncbi_analysis.md) — detailed methods report
- [docs/coexpression_network_methods.md](docs/coexpression_network_methods.md) — co-expression network & guilt-by-association methods

### Reference (`data/reference/`)
- [ref_gene_table.csv](data/reference/ref_gene_table.csv) — gene_id, chrom, coords, biotype, product
- [tx2gene.tsv](data/reference/tx2gene.tsv) — transcript→gene map
- [fgsg_go_map.tsv](data/reference/fgsg_go_map.tsv) — gene→GO (from blast2go GAF)
- [go_term_names.tsv](data/reference/go_term_names.tsv) — GO ID→name (EBI QuickGO)

### Expression matrices (`data/matrices/`) — 13,725 gene models × 12 libraries
- [raw_counts.csv](data/matrices/raw_counts.csv) — featureCounts raw counts
- [tpm_matrix.csv](data/matrices/tpm_matrix.csv) — StringTie TPM
- [fpkm_matrix.csv](data/matrices/fpkm_matrix.csv) — StringTie FPKM
- [normalized_counts.csv](data/matrices/normalized_counts.csv) — DESeq2 median-of-ratios
- [vst_matrix.csv](data/matrices/vst_matrix.csv) — variance-stabilized (heatmap input)

### QC (`results/qc/`)
- [fastp_summary.csv](results/qc/fastp_summary.csv) · [hisat2_alignment_summary.csv](results/qc/hisat2_alignment_summary.csv)

### Differential expression
- DESeq2 (raw counts): [`results/de_deseq2/`](results/de_deseq2/) — 6 per-contrast tables + [DE_summary_counts.csv](results/de_deseq2/DE_summary_counts.csv) + [all_significant_DEGs_long.csv](results/de_deseq2/all_significant_DEGs_long.csv)
- limma-trend (TPM): [`results/de_tpm/`](results/de_tpm/) — 6 tables + [DE_summary_tpm.csv](results/de_tpm/DE_summary_tpm.csv)
- limma-trend (FPKM): [`results/de_fpkm/`](results/de_fpkm/) — 6 tables + [DE_summary_fpkm.csv](results/de_fpkm/DE_summary_fpkm.csv)
- Cross-method: [results/concordance/method_concordance.csv](results/concordance/method_concordance.csv)

### GO enrichment (`results/go/`)
- [GO_enrichment_all_contrasts.csv](results/go/GO_enrichment_all_contrasts.csv) · [GO_enrichment_significant.csv](results/go/GO_enrichment_significant.csv) + 6 per-contrast tables

### Genes of interest (`results/genes_of_interest/`)
- [GOI_expression_DESeq2_6contrasts.xlsx](results/genes_of_interest/GOI_expression_DESeq2_6contrasts.xlsx) — 6 contrasts, DESeq2
- [GOI_expression_template_3contrasts_FPKM.xlsx](results/genes_of_interest/GOI_expression_template_3contrasts_FPKM.xlsx) — CLC-style template, 3 contrasts, FPKM
- [GOI_expression_all_methods_6contrasts.xlsx](results/genes_of_interest/GOI_expression_all_methods_6contrasts.xlsx) — 6 contrasts × 3 methods
- [GOI_expression_long_all_methods.csv](results/genes_of_interest/GOI_expression_long_all_methods.csv) · [GOI_overlap_by_category.csv](results/genes_of_interest/GOI_overlap_by_category.csv)

### Venn membership (`results/venn/`)
- [venn_treatment_vs_control_membership.csv](results/venn/venn_treatment_vs_control_membership.csv)

### Co-expression network (`results/coexpression_network/`)
- [gene_modules.csv](results/coexpression_network/gene_modules.csv) — WGCNA module assignment for all 12,723 genes
- [coexpression_edges_goi_anchored.csv](results/coexpression_network/coexpression_edges_goi_anchored.csv) — 32,730 GOI→candidate edges (TOM ≥ 0.40)
- [candidate_genes.csv](results/coexpression_network/candidate_genes.csv) — 2,982 candidates with partner counts, TOM, module, product, GO
- [candidate_ranking_combined.csv](results/coexpression_network/candidate_ranking_combined.csv) — combined cross-treatment ranking
- [candidate_ranking_black.csv](results/coexpression_network/candidate_ranking_black.csv) · [candidate_ranking_dark.csv](results/coexpression_network/candidate_ranking_dark.csv) · [candidate_ranking_low.csv](results/coexpression_network/candidate_ranking_low.csv) — per-treatment GBA rankings
- [goi_treatment_membership.csv](results/coexpression_network/goi_treatment_membership.csv) — GOI DEG membership per treatment
- [goi_category_map.csv](results/coexpression_network/goi_category_map.csv) — gene→functional category for all 257 GOI
- [candidate_top40_GO_enrichment.csv](results/coexpression_network/candidate_top40_GO_enrichment.csv) — GO over-representation of the top-40 candidates (13 terms; hits, fold-enrichment, padj, candidate FGSG IDs per term)
- [candidate_GO_subnetwork_membership.csv](results/coexpression_network/candidate_GO_subnetwork_membership.csv) — per-GO-term subgraph membership (gene, role: GOI/candidate_seed/candidate_other, category)

### Figures (`figures/`)
All 18 PNGs referenced above.

## Tool versions

fastp 1.1.0 · HISAT2 2.2.2 · samtools 1.22.1 · subread/featureCounts 2.1.1 ·
StringTie 3.0.3 · DESeq2 (R 4.5.3) + apeglm · limma · pheatmap

## License

Released under the MIT License — see [LICENSE](LICENSE).
