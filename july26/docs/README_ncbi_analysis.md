# Fusarium graminearum PH-1 — Light-Condition RNA-seq Re-analysis (NCBI/GenBank reference)

Differential expression of *Fusarium graminearum* PH-1 grown under four light regimes,
re-run against the **NCBI RefSeq genome and annotation** (assembly GCF_000240135.3,
ASM24013v3) so that all results carry native **FGSG_ gene IDs**. This supersedes the
earlier funannotate-based build; the pipeline here quantifies expression three ways
(raw counts, TPM, FPKM), runs differential expression on all three, clusters the raw
count-derived expression, and performs GO enrichment against the blast2go annotation.

## Design

| Condition | Code | Samples | Description |
|-----------|------|---------|-------------|
| Control   | H | H13_S7, H14_S8, H15_S9    | half/half — 12 h light / 12 h dark (reference level) |
| Black light | B | B10_S1, B11_S2, B12_S3  | black light |
| Dark light  | D | D10_S4, D11_S5, D12_S6  | dark light |
| Low light   | L | L13_S10, L14_S11, L15_S12 | low light |

3 biological replicates per condition, 12 paired-end libraries (~101 bp).
All **six pairwise contrasts** were tested.

## Pipeline

1. **Reference** — HISAT2 index built from the NCBI genome (19 sequences: 4 chromosomes
   NC_026474–77 + 15 unplaced contigs). Annotation: 13,725 genes (13,312 protein-coding,
   314 tRNA, 88 rRNA, plus pseudogenes), all FGSG_-prefixed.
2. **Trimming** — fastp (adapter auto-detect, Q20, min length 36). ~103.7 M read pairs
   retained (>98.5% pass, Q30 ≈ 93%).
3. **Alignment** — HISAT2 (`--rna-strandness RF --dta`), sorted/indexed with samtools.
   Overall alignment 97.3–97.9%. *B12_S3 shows 30% concordant multi-mapping — a known
   property of that library, retained; its multimappers are excluded from counting.*
4. **Quantification**
   - **Raw counts** — featureCounts (`-p --countReadPairs -s 2`, gene level, exon features).
   - **TPM & FPKM** — StringTie (`-e -B --rf`) gene abundances.
   - Output: 13,725 genes × 12 samples for each; gene sets identical across the three matrices;
     TPM columns sum to exactly 1e6.
5. **QC** — DESeq2 VST → PCA (PC1 = 80%, PC2 = 10%), Spearman sample correlation
   (r ≈ 0.99 within condition, ≈ 0.93 between), library size and gene-detection profiles.
   Samples separate cleanly by condition; replicates cluster tightly.
6. **Differential expression** — three parallel analyses, threshold **padj < 0.05 & |log2FC| ≥ 1**:
   - **Raw counts → DESeq2** (Wald test, apeglm LFC shrinkage).
   - **TPM → limma-trend** on log2(TPM+1) (robust empirical-Bayes, trend=TRUE).
   - **FPKM → limma-trend** on log2(FPKM+1).
7. **Heatmaps** — VST-transformed expression (derived from raw counts), row z-scored,
   Ward hierarchical clustering on both axes.
8. **GO enrichment** — hypergeometric over-representation against the blast2go GAF
   (`blast2go_FgraPH1_FGSG.gaf`, 8,961 annotated genes, 36,729 gene–GO pairs),
   BP/MF/CC tested separately, Benjamini–Hochberg FDR, background = genes tested in each
   contrast, requiring ≥2 DEGs per term. GO term names resolved from EBI QuickGO.

## Differential expression — DEG counts (padj<0.05, |log2FC|≥1)

**DESeq2 (raw counts)**

| Contrast | tested | sig | up in first | up in second |
|----------|-------:|----:|------------:|-------------:|
| dark vs black | 11237 | 1989 | 1155 | 834 |
| control vs black | 10990 | 753 | 426 | 327 |
| low vs black | 11483 | 2721 | 1503 | 1218 |
| control vs dark | 10990 | 1219 | 572 | 647 |
| low vs dark | 10990 | 631 | 254 | 377 |
| low vs control | 11237 | 1639 | 853 | 786 |

3,525 unique DEGs across all contrasts; 479 shared across the three light-vs-black
contrasts; 42 DE in all six contrasts.

**limma-trend (TPM)** — DEG totals: 1965 / 780 / 2877 / 1116 / 617 / 1634 (same contrast order).
**limma-trend (FPKM)** — DEG totals: 1743 / 639 / 2681 / 950 / 526 / 1436.

### Cross-method concordance

The three quantifications give closely matched DEG counts and highly correlated
fold-changes (pooled Pearson r = 0.94 between DESeq2 and limma-TPM log2FC).
Three-way significant-gene overlap ranges 48–65% per contrast (Jaccard-style; the
union/intersection breakdown is in `method_concordance.csv`). The rank order of
contrasts by DEG count is identical across methods (low-vs-black largest,
low-vs-dark smallest).

## GO enrichment

168 GO terms significant at padj < 0.05 across the six contrasts. The light-response
signature is dominated by **oxidation-reduction / oxidoreductase activity**
(cytochrome P450 / monooxygenase, heme & iron binding), **oxidative-stress detoxification**
(catalase, catalase-peroxidase, hydrogen-peroxide catabolism, thioredoxin/TSA
antioxidant enzymes), and **amino-acid / small-molecule transport** (amino-acid and
lactose permeases). This is biologically consistent with light-driven redox and
photo-oxidative-stress metabolism in *Fusarium*.

*Note on obsolete terms:* the 2023 blast2go annotation uses `GO:0055114`
(oxidation-reduction process), which has since been made obsolete in the current Gene
Ontology; it is reported as-is because it is the annotation actually present in the GAF.

## Files

### Expression matrices (`matrices/`)
- `raw_counts.csv` — featureCounts raw integer counts (gene_id, biotype, product, 12 samples)
- `tpm_matrix.csv`, `fpkm_matrix.csv` — StringTie TPM / FPKM
- `vst_matrix.csv` — DESeq2 variance-stabilized expression (used for QC & heatmaps)
- `normalized_counts.csv` — DESeq2 median-of-ratios normalized counts

### Differential expression
- `de_counts/DE_<contrast>.csv` — DESeq2, 6 contrasts (+ `DE_summary_counts.csv`, `all_significant_DEGs_long.csv`)
- `de_tpm/`, `de_fpkm/` — limma-trend equivalents
- `method_concordance.csv` — cross-method overlap per contrast

### QC & figures (`figures/`)
- `qc_alignment_rates.png`, `qc_pca.png`, `qc_correlation.png`, `qc_libsize_detection.png`
- `volcano_counts_grid.png`, `ma_counts_grid.png` — DESeq2 volcano & MA (6 contrasts)
- `method_concordance.png` — DEG counts & log2FC agreement across methods
- `heatmap_all_genes_clustered.png` — all 12,723 expressed gene models
- `heatmap_DEGs_clustered.png` — 3,525 DEG union
- `heatmap_top50_DEGs_labeled.png` — top-50 variable DEGs with FGSG labels
- `go_dotplot.png` — top enriched GO terms per contrast

### GO (`go/`)
- `GO_enrichment_all_contrasts.csv`, `GO_enrichment_significant.csv`, `GO_<contrast>.csv`

### Reference (`ref/`)
- `ref_gene_table.csv`, `tx2gene.tsv`, `fgsg_go_map.tsv`, `go_term_names.tsv`

## Tool versions

fastp 1.1.0 · HISAT2 2.2.2 · samtools 1.22.1 · subread/featureCounts 2.1.1 ·
StringTie 3.0.3 · DESeq2 (R 4.5.3) + apeglm · limma (bioconductor) · GO names via EBI QuickGO.
