# *Fusarium graminearum* PH-1 RNA-seq — Four Light Conditions

RNA-seq processing and differential expression for *Fusarium graminearum*
PH-1 grown under **four light conditions** — black, dark, high, and low
light — with **3 biological replicates each (12 libraries)** and **all six
pairwise differential-expression contrasts**.

**Pipeline:** Salmon (selective alignment → TPM + counts) + HISAT2 (spliced
genome alignment → BAMs) + StringTie (FPKM/TPM) → DESeq2 (pairwise DE).
Runs natively on Apple Silicon (osx-arm64).

## Repository layout

| Folder | Contents |
|---|---|
| [`figures/`](figures) | All QC, DE, expression-heatmap, and Venn figures (PNG) |
| [`matrices/`](matrices) | Expression matrices, normalized counts, VST, PCA, correlation, mapping summary, FGSG cross-reference |
| [`de_results/`](de_results) | Per-contrast DE tables, combined long table, DE summary |
| [`scripts/`](scripts) | DESeq2 pairwise DE + QC scripts, sample-sheet template |

All gene-level tables (expression matrices and DE results) now carry an
`FGSG_id` column immediately after `gene_id`, mapping each funannotate
`FGRAM000v1_` model to its legacy Broad `FGSG_` identifier — see
[FGSG cross-reference](#fgsg-cross-reference) below.

---

## Experimental design

| Condition | Replicates |
|---|---|
| black_light | B10_S1, B11_S2, B12_S3 |
| dark_light  | D10_S4, D11_S5, D12_S6 |
| high_light  | H13_S7, H14_S8, H15_S9 |
| low_light   | L13_S10, L14_S11, L15_S12 |

- **Genome:** *F. graminearum* PH-1, assembly ASM24013v3 (19 sequences: chr1–4 + 15 contigs)
- **Annotation:** funannotate GFF3 — 13,377 genes / 13,569 mRNA / 301 tRNA
- **Reads:** paired-end, ~100 bp, 6–10M pairs per library

## Tool versions

salmon 2.3.4 · hisat2 2.2.2 · samtools 1.22.1 · fastp 1.1.0 ·
stringtie 3.0.3 · gffread 0.12.9 · subread/featureCounts 2.1.1 ·
DESeq2 1.50.2 · tximport 1.38.2 · apeglm (R 4.5.3)

## Workflow

1. **Reference prep** — `gffread` extracted 13,870 transcripts; tx2gene map (13,870 tx → 13,377 genes); Salmon decoy-aware index (whole genome as 19 decoys, k=31); HISAT2 genome index + splice sites.
2. **Trimming** — `fastp --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 36`. >98% reads passed, Q30 ~93%, GC ~51% across all 12.
3. **Salmon quant** — `-l A --seqBias --gcBias -g tx2gene`. Library type auto-detected **ISR** (stranded, fr-firststrand).
4. **HISAT2 align** — `--rna-strandness RF --dta` + known splice sites → sorted/indexed BAMs (for IGV).
5. **StringTie** — `-e -G Fgram.gtf --rf` → per-gene FPKM & TPM.
6. **Matrices** — Salmon gene counts (DESeq2 input via tximport) + TPM; StringTie FPKM/TPM.
7. **QC** — library size, gene detection, replicate correlation, VST PCA.
8. **Pairwise DE** — [`scripts/run_pairwise_DE.R`](scripts/run_pairwise_DE.R): all 6 DESeq2 Wald contrasts, apeglm LFC shrinkage, padj < 0.05 & |log2FC| ≥ 1.

---

## Quality control

### Sample structure (VST PCA)

![VST PCA of all 12 libraries]({{artifact:art_d9cb2822-6a13-4259-a6c2-812152883ddc}})

The four conditions separate cleanly with tight replicate grouping. PC1 = 81%,
PC2 = 11% of variance (top 500 variable genes). Within-condition replicate
correlation r ≈ 0.997 vs r ≈ 0.955 between conditions.

### Replicate correlation

![Sample correlation heatmap]({{artifact:art_a073a21d-a0ef-45dd-977c-cd3e72b0b372}})

### Library size & gene detection

![Library size and gene detection]({{artifact:art_3c0c496f-4ab1-4737-864c-8b70b58fa93c}})

9,100–9,600 genes detected (TPM > 1) per library, uniform across conditions.

### Mapping summary

Full table: [`matrices/mapping_alignment_summary_12samples.csv`](matrices/mapping_alignment_summary_12samples.csv)

| sample | Salmon %mapped | HISAT2 overall | HISAT2 multimap |
|---|---|---|---|
| B10_S1  | 80.6% | 97.5% | 2.2% |
| B11_S2  | 82.6% | 97.5% | 0.5% |
| B12_S3  | 55.2% | 97.4% | 30.0% |
| D10_S4  | 82.6% | 97.7% | 1.2% |
| D11_S5  | 83.3% | 97.9% | 1.8% |
| D12_S6  | 81.9% | 97.8% | 1.6% |
| H13_S7  | 81.1% | 97.5% | 0.7% |
| H14_S8  | 83.1% | 97.6% | 0.6% |
| H15_S9  | 81.2% | 97.7% | 1.5% |
| L13_S10 | 85.2% | 97.7% | 0.5% |
| L14_S11 | 85.3% | 97.7% | 0.5% |
| L15_S12 | 84.1% | 97.8% | 1.2% |

> **Note on B12_S3:** aligns to the genome normally (97%) but has a low Salmon
> transcriptome rate (55%) and 30% multimapping — indicating more
> rRNA/repetitive/intergenic content. Despite this it clusters tightly with
> B10/B11 on the PCA (expression profile intact), so it is retained.

---

## Differential expression — 6 pairwise contrasts

DESeq2 Wald test, apeglm-shrunk LFC. Significance: **padj < 0.05 AND
|log2FC| ≥ 1**. In each `DE_<B>_vs_<A>.csv`, positive log2FoldChange = UP in
the first-named (B) condition. 9,693 genes tested per contrast.

![DEG counts across all six contrasts]({{artifact:art_17736282-e66d-437c-8ef3-8874ee37f87a}})

| Contrast (B vs A) | Sig DEGs | Up in B | Down in B | Table |
|---|---|---|---|---|
| dark vs black  | 2,020 | 1,195 | 825   | [csv](de_results/DE_dark_light_vs_black_light.csv) |
| high vs black  |   857 |   518 | 339   | [csv](de_results/DE_high_light_vs_black_light.csv) |
| low vs black   | 2,590 | 1,475 | 1,115 | [csv](de_results/DE_low_light_vs_black_light.csv) |
| high vs dark   | 1,308 |   585 | 723   | [csv](de_results/DE_high_light_vs_dark_light.csv) |
| low vs dark    |   698 |   298 | 400   | [csv](de_results/DE_low_light_vs_dark_light.csv) |
| low vs high    | 1,703 |   910 | 793   | [csv](de_results/DE_low_light_vs_high_light.csv) |

Summary table: [`de_results/DE_summary.csv`](de_results/DE_summary.csv) ·
Combined long table (every significant gene×contrast row):
[`de_results/all_significant_DEGs_long.csv`](de_results/all_significant_DEGs_long.csv)

### Global patterns

- **3,330 unique genes** are DE in ≥1 contrast.
- **561 core light-response genes** are shared across all three light-vs-black contrasts (dark/high/low each vs black) — a candidate core transcriptional response to illumination.
- **43 genes** are DE in all six contrasts (distinguish every condition pair).
- Low light drives the largest response vs black (2,590 DEGs); high light the smallest (857). Low-vs-dark is the most similar pair (698 DEGs).

### Volcano plots

| | |
|---|---|
| ![dark vs black]({{artifact:art_3d33ce01-d848-40b2-a14c-c6a6846db3ab}}) | ![high vs black]({{artifact:art_e7ffb32a-8234-4a44-9f1d-7f33f6cd5053}}) |
| ![low vs black]({{artifact:art_f564d7ed-2e29-48db-bbca-54ed2e76b605}}) | ![high vs dark]({{artifact:art_6b50752d-7bed-496d-93f8-b89cb1a5ea2a}}) |
| ![low vs dark]({{artifact:art_1165077c-ebce-48e7-bfdd-374b51e65ee4}}) | ![low vs high]({{artifact:art_b1f48d56-8035-4818-9734-3144639265a5}}) |

### MA plots

| | |
|---|---|
| ![dark vs black]({{artifact:art_aae1ba2d-92e8-4d08-9bfa-743b411db44d}}) | ![high vs black]({{artifact:art_54fb6c5d-6e98-4ef9-aa14-b7f5e46dc26f}}) |
| ![low vs black]({{artifact:art_cfa79b9e-bbdf-44f9-a05b-5e2458aa4fec}}) | ![high vs dark]({{artifact:art_64698652-6380-46c7-b3d6-0ae336f54a0d}}) |
| ![low vs dark]({{artifact:art_1c869941-8ea5-40f7-8999-043b3f168e10}}) | ![low vs high]({{artifact:art_58256273-2fea-4016-8b09-6f182ebabfd6}}) |

---

## Expression heatmaps

Replicates merged by **per-condition mean TPM**; each gene row is a
**z-score of log₂(mean TPM + 1)** across the four conditions, hierarchically
clustered (average linkage, Euclidean distance). Blue = below a gene's mean,
red = above. Columns ordered black · dark · high · low.

### Top 60 most variable DEGs (labeled)

![Top 60 variable DEGs]({{artifact:art_6d2cd4d2-b212-4206-8e73-273584f5493c}})

The 60 DEGs with the largest cross-condition variance, row-labeled with
`FGSG_id` where available. Condition-specific up-regulated blocks are clearly
resolved (e.g. a large low-light-specific cluster).

### All 3,330 DEGs (global structure)

![All DEGs]({{artifact:art_1fd0e589-3658-43e8-ac5e-e5f338e00018}})

Every gene DE in ≥1 contrast (row labels omitted for density). Confirms the
four conditions partition into coherent co-expression blocks.

---

## Condition overlap — Venn diagrams

Genes are called **"expressed" in a condition when the mean expression across
its 3 replicates meets a cutoff**, then the four resulting gene sets are
intersected. This is an expression-occupancy view (where transcription is
present), distinct from the DE analysis above (where transcription *changes*).
Four cutoff instances are provided:

| | |
|---|---|
| ![Venn TPM≥1]({{artifact:art_4c12296c-a4b5-48d1-a383-ade24943e65b}}) | ![Venn TPM≥5]({{artifact:art_ee09bf19-9d6f-44a6-8188-49cb83d397a1}}) |
| ![Venn TPM≥10]({{artifact:art_33395130-a0bb-452a-8ed2-107e93f5425d}}) | ![Venn FPKM≥5]({{artifact:art_8871363d-e194-4176-8f44-7fabd728b363}}) |

The large majority of genes are expressed in **all four** conditions
(core transcriptome), with the biggest condition-restricted and pairwise
compartments consistently involving **low light** and the **dark–low** pair.
Region counts for all 15 subsets × 4 cutoffs:
[`matrices/venn_region_counts.csv`](matrices/venn_region_counts.csv).

| Region | TPM≥1 | TPM≥5 | TPM≥10 | FPKM≥5 |
|---|--:|--:|--:|--:|
| **All four (core)** | 8,811 | 7,051 | 5,878 | 6,453 |
| Dark & Low only | 212 | 310 | 359 | 325 |
| Dark/High/Low | 252 | 315 | 326 | 274 |
| Low only | 156 | 259 | 319 | 283 |
| Black & High only | 99 | 92 | 70 | 85 |
| Total expressed (≥1 condition) | 9,999 | 8,435 | 7,402 | 7,874 |

---

## FGSG cross-reference

The funannotate annotation uses custom `FGRAM000v1_` gene IDs. Legacy Broad
Institute `FGSG_` identifiers were parsed **authoritatively from the GFF3
`product=` field** (the only species-specific source of this mapping) and
folded back into every gene-level table as an `FGSG_id` column.

Full mapping table:
[`matrices/fgram_fgsg_crossref.csv`](matrices/fgram_fgsg_crossref.csv)
— columns `gene_id, FGSG_id, product, fgsg_shared_with_other_gene`.

- **1,903 of 13,377 genes** carry an `FGSG_id`; the remaining genes have no
  legacy equivalent in this annotation (blank in the column). No external
  extension is possible because the `FGRAM000v1_` scheme is custom to this build.
- The mapping is **1:1** for all matched genes — **0 conflicts**.
- **5 FGSG IDs are each shared by 2 consecutive FGRAM genes** (legacy models
  that funannotate split into two). Both genes are kept and flagged
  `fgsg_shared_with_other_gene = True`.

---

## Data files

All gene-level matrices below carry an `FGSG_id` column after `gene_id`.

- [`gene_counts_salmon.csv`](matrices/gene_counts_salmon.csv) — gene × 12 samples estimated counts (DESeq2 input)
- [`gene_TPM_salmon.csv`](matrices/gene_TPM_salmon.csv) — gene × 12 samples TPM (Salmon)
- [`stringtie_TPM_matrix.csv`](matrices/stringtie_TPM_matrix.csv), [`stringtie_FPKM_matrix.csv`](matrices/stringtie_FPKM_matrix.csv) — with gene coordinates
- [`normalized_counts_allsamples.csv`](matrices/normalized_counts_allsamples.csv) — DESeq2 median-of-ratios normalized
- [`vst_allsamples.csv`](matrices/vst_allsamples.csv) — variance-stabilized values
- [`pca_allsamples.csv`](matrices/pca_allsamples.csv), [`sample_correlation_allsamples.csv`](matrices/sample_correlation_allsamples.csv)
- [`fgram_fgsg_crossref.csv`](matrices/fgram_fgsg_crossref.csv) — FGRAM000v1_ ↔ FGSG_ mapping + product + split-gene flag
- [`venn_region_counts.csv`](matrices/venn_region_counts.csv) — Venn region counts (15 subsets × 4 cutoffs)

### Scripts ([`scripts/`](scripts))
- [`run_pairwise_DE.R`](scripts/run_pairwise_DE.R) — automated all-pairwise DESeq2 (Wald + apeglm), volcano + MA per contrast
- [`qc_allsamples.R`](scripts/qc_allsamples.R) — VST, PCA, correlation, library QC
- [`condition_metadata_template.tsv`](scripts/condition_metadata_template.tsv) — sample-sheet template

Re-run all DE:
```bash
Rscript scripts/run_pairwise_DE.R samples.tsv salmon/ tx2gene.tsv de/
```

## Reproducibility

- conda env `rnaseq` (osx-arm64): alignment/quant tools
- conda env `deseq` (R 4.5.3): DESeq2 stack
- Sorted BAMs (`<sample>.sorted.bam`, for IGV) are produced by the pipeline but not included here due to size.
