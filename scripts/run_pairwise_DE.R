#!/usr/bin/env Rscript
# =====================================================================
# Pairwise differential expression across Fusarium light conditions
# Salmon -> tximport -> DESeq2, all pairwise Wald contrasts
#
# Usage:
#   Rscript run_pairwise_DE.R <samples.tsv> <salmon_dir> <tx2gene.tsv> <out_dir>
#
# samples.tsv columns (tab-sep, header):
#   sample  condition  R1  R2
# Add rows for dark_light / high_light / low_light and rerun to get the
# full pairwise matrix automatically.
# =====================================================================
suppressMessages({
  library(DESeq2); library(tximport); library(ggplot2); library(ggrepel)
})

args <- commandArgs(trailingOnly=TRUE)
samples_tsv <- ifelse(length(args)>=1, args[1],
  "/Users/carsonandorf/Documents/ClaudeScience/expression_ph1/analysis/ref/samples.tsv")
salmon_dir  <- ifelse(length(args)>=2, args[2],
  "/Users/carsonandorf/Documents/ClaudeScience/expression_ph1/analysis/salmon")
tx2gene_tsv <- ifelse(length(args)>=3, args[3],
  "/Users/carsonandorf/Documents/ClaudeScience/expression_ph1/analysis/ref/tx2gene.tsv")
out_dir     <- ifelse(length(args)>=4, args[4],
  "/Users/carsonandorf/Documents/ClaudeScience/expression_ph1/analysis/de")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# thresholds
ALPHA <- 0.05      # adjusted p-value cutoff
LFC   <- 1.0       # |log2FC| cutoff for "significant" calls (>=2-fold)

# ---- load sample sheet ----
ss <- read.table(samples_tsv, sep="\t", header=TRUE, stringsAsFactors=FALSE)
ss$condition <- factor(ss$condition)
rownames(ss) <- ss$sample
conds <- levels(ss$condition)
cat("Conditions:", paste(conds, collapse=", "),
    "| samples:", nrow(ss), "\n")

# ---- tximport ----
files <- file.path(salmon_dir, ss$sample, "quant.sf")
names(files) <- ss$sample
stopifnot(all(file.exists(files)))
tx2gene <- read.table(tx2gene_tsv, sep="\t",
                      col.names=c("tx","gene"), stringsAsFactors=FALSE)
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

# ---- DESeq dataset ----
# With a single condition, DESeq2 requires design=~1 (no contrasts possible).
single_cond <- length(conds) < 2
design_formula <- if (single_cond) ~1 else ~condition
dds <- DESeqDataSetFromTximport(txi, colData=ss, design=design_formula)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep,]
cat("Genes after prefilter (>=10 counts in >=2 samples):", nrow(dds), "\n")

if (single_cond) {
  # Validate normalization/dispersion end-to-end, export normalized counts, then stop.
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, fitType="local", quiet=TRUE)
  norm_counts <- counts(dds, normalized=TRUE)
  write.csv(round(norm_counts,2), file.path(out_dir,"normalized_counts.csv"))
  cat("\n[Only one condition present:", conds[1], "] DE contrasts require >=2 conditions.\n")
  cat("Normalization + dispersion validated end-to-end (", nrow(dds), "genes ).\n")
  cat("Add dark_light / high_light / low_light rows to samples.tsv and rerun\n")
  cat("to get all pairwise contrasts automatically.\n")
  quit(save="no", status=0)
}

# fitType="local" is robust for small n; parametric can fail
dds <- DESeq(dds, fitType="local", quiet=TRUE)

# Export normalized counts
norm_counts <- counts(dds, normalized=TRUE)
write.csv(round(norm_counts,2), file.path(out_dir,"normalized_counts.csv"))

pairs <- combn(conds, 2, simplify=FALSE)
summary_rows <- list()
for (pr in pairs) {
  a <- pr[1]; b <- pr[2]              # contrast: b vs a (a = reference)
  tag <- paste0(b, "_vs_", a)
  res <- results(dds, contrast=c("condition", b, a), alpha=ALPHA)
  # shrink LFC for ranking/plots (apeglm needs coef; use ashr-free normal via lfcShrink type="normal")
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[order(res_df$padj),
                   c("gene_id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
  write.csv(res_df, file.path(out_dir, paste0("DE_", tag, ".csv")), row.names=FALSE)

  sig <- subset(res_df, !is.na(padj) & padj < ALPHA & abs(log2FoldChange) >= LFC)
  up  <- sum(sig$log2FoldChange > 0); dn <- sum(sig$log2FoldChange < 0)
  summary_rows[[tag]] <- data.frame(contrast=tag, tested=nrow(res_df),
                                    sig=nrow(sig), up_in_b=up, down_in_b=dn)

  # ---- volcano ----
  v <- res_df[!is.na(res_df$padj),]
  v$sig <- with(v, ifelse(padj<ALPHA & abs(log2FoldChange)>=LFC,
                          ifelse(log2FoldChange>0,"up","down"),"ns"))
  p <- ggplot(v, aes(log2FoldChange, -log10(padj), color=sig)) +
    geom_point(size=0.7, alpha=0.6) +
    scale_color_manual(values=c(up="#b2182b", down="#2166ac", ns="grey75")) +
    geom_vline(xintercept=c(-LFC,LFC), linetype="dashed", linewidth=0.3) +
    geom_hline(yintercept=-log10(ALPHA), linetype="dashed", linewidth=0.3) +
    labs(title=paste0(b, " vs ", a),
         subtitle=paste0(nrow(sig)," DE genes (padj<",ALPHA,", |LFC|>=",LFC,")"),
         x="log2 fold change", y="-log10 adjusted p") +
    theme_bw(base_size=9) + theme(legend.title=element_blank())
  ggsave(file.path(out_dir, paste0("volcano_", tag, ".png")),
         plot=p, width=4.5, height=3.8, dpi=200)

  # ---- MA ----
  png(file.path(out_dir, paste0("MA_", tag, ".png")), width=1000, height=800, res=200)
  plotMA(res, main=paste0(b," vs ",a), alpha=ALPHA); dev.off()
}

de_summary <- do.call(rbind, summary_rows)
write.csv(de_summary, file.path(out_dir,"DE_summary.csv"), row.names=FALSE)
cat("\n=== DE summary ===\n"); print(de_summary, row.names=FALSE)
cat("\nWrote per-contrast tables, volcano + MA plots to:", out_dir, "\n")
