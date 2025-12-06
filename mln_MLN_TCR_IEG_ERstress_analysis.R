# ============================================================
# Part 2 — MLN-only IEG / Nur77 / ER-stress analysis
# Dataset: mLN iNKT cells (WT_PBS, WT_aGC, KO_PBS)
# Author: Pumla Bhekiwe Manyatsi
#
# This script is a continuation of the MLN vs Spleen DE analysis:
#   - Part 1: MLN vs Spleen (mln_spleen_DE_analysis.Rmd)
#   - Part 2 (this script): within-MLN comparison of
#       WT + PBS vs WT + αGalCer vs CD1d-KO + PBS
# ============================================================

# -------------------------------
# STEP 0 — Setup
# -------------------------------

rm(list = ls(all.names = TRUE))

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(pheatmap)
  library(biomaRt)
  library(dplyr)
  library(tidyr)
})

# For GitHub, prefer relative paths instead of /archive/...
# Put your STAR outputs in: data/MLN_ReadsPerGene/
base_dir <- "data/MLN_ReadsPerGene"
out_dir  <- "results/MLN_MLN_TCR_IEG_ERstress"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------------------------------
# STEP 1 — Define MLN-only groups
# -------------------------------

# WT MLN + PBS (baseline)
wt_pbs_files <- c(
  "SRR25987863_pass_ReadsPerGene.out.tab",  # GSM7770022
  "SRR25987867_pass_ReadsPerGene.out.tab",  # GSM7770021
  "SRR25987871_pass_ReadsPerGene.out.tab"   # GSM7770020
)

# WT MLN + αGalCer (TCR-stimulated)
wt_agc_files <- c(
  "SRR25987840_pass_ReadsPerGene.out.tab",  # GSM7770027
  "SRR25987836_pass_ReadsPerGene.out.tab"   # GSM7770028
)

# CD1d KO MLN + PBS (negative control)
ko_pbs_files <- c(
  "SRR25987820_ReadsPerGene.out.tab",  # GSM7770032
  "SRR25987816_ReadsPerGene.out.tab",  # GSM7770033
  "SRR25987812_ReadsPerGene.out.tab"   # GSM7770034
)

# OPTIONAL: CD1d KO MLN + αGalCer (not used yet)
ko_agc_files <- c(
  "SRR25987892_pass_ReadsPerGene.out.tab",  # GSM7770038
  "SRR25987888_pass_ReadsPerGene.out.tab",  # GSM7770039
  "SRR25987884_pass_ReadsPerGene.out.tab"   # GSM7770040
)

# Combined list (WT_PBS + WT_aGC + KO_PBS)
target_files <- c(wt_pbs_files, wt_agc_files, ko_pbs_files)
file_list_mln <- file.path(base_dir, target_files)

cat("MLN file list:\n")
print(file_list_mln)
cat("Files exist?\n")
print(file.exists(file_list_mln))

# -------------------------------
# STEP 2 — Count matrix + metadata
# -------------------------------

# Keep only essential objects
rm(list = setdiff(ls(), c("file_list_mln", "target_files", "base_dir", "out_dir")))

read_counts_mln <- function(file) {
  df <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  df <- df[grep("^ENS", df$V1), c(1, 2)]  # Ensembl only, column 2 = unstranded
  sample_name <- gsub("(_pass)?_ReadsPerGene\\.out\\.tab$", "", basename(file))
  colnames(df) <- c("GeneID", sample_name)
  df
}

data_list_mln <- lapply(file_list_mln, read_counts_mln)

combined_counts_mln <- Reduce(
  function(x, y) merge(x, y, by = "GeneID", all = TRUE),
  data_list_mln
)
combined_counts_mln[is.na(combined_counts_mln)] <- 0

count_data_mln <- as.matrix(combined_counts_mln[, -1])
rownames(count_data_mln) <- combined_counts_mln$GeneID

sample_names_mln <- colnames(count_data_mln)

## file_list_mln order:
##  1–3: WT + PBS
##  4–5: WT + αGalCer
##  6–8: KO + PBS

Genotype <- factor(c(
  rep("WT", 5),
  rep("KO", 3)
))

Treatment <- factor(c(
  rep("PBS", 3),  # WT PBS
  rep("aGC", 2),  # WT αGalCer
  rep("PBS", 3)   # KO PBS
))

Group <- interaction(Genotype, Treatment, sep = "_")  # WT_PBS, WT_aGC, KO_PBS

sample_info_mln <- data.frame(
  Sample    = sample_names_mln,
  Genotype  = Genotype,
  Treatment = Treatment,
  Group     = Group,
  row.names = sample_names_mln
)

print(sample_info_mln)

# -------------------------------
# STEP 3 — edgeR + PCA
# -------------------------------

dge <- DGEList(counts = count_data_mln, samples = sample_info_mln)
keep <- filterByExpr(dge, group = dge$samples$Group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

logCPM <- cpm(dge, log = TRUE, prior.count = 1)

write.csv(logCPM, file.path(out_dir, "MLN_logCPM_matrix.csv"))
write.csv(dge$samples, file.path(out_dir, "MLN_sample_metadata.csv"))

pca <- prcomp(t(logCPM), scale. = TRUE)
pca_df <- data.frame(
  PC1    = pca$x[,1],
  PC2    = pca$x[,2],
  Group  = dge$samples$Group,
  Sample = rownames(dge$samples)
)

p <- ggplot(pca_df, aes(PC1, PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.6, size = 3) +
  theme_bw() +
  ggtitle("PCA: mLN (WT_PBS vs WT_aGC vs KO_PBS)") +
  theme(plot.title = element_text(hjust = 0.5))

print(p)
ggsave(file.path(out_dir, "PCA_MLN_TCR_groups.pdf"), p, width = 6, height = 5)

# -------------------------------
# STEP 4 — Ensembl → symbol
# -------------------------------

ens_ids <- sub("\\..*", "", rownames(logCPM))

mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

annot_tbl <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters    = "ensembl_gene_id",
  values     = unique(ens_ids),
  mart       = mart
)

anno_df <- data.frame(
  ensembl_gene_id = ens_ids,
  stringsAsFactors = FALSE
) %>%
  left_join(annot_tbl, by = "ensembl_gene_id")

gene_symbols <- ifelse(
  is.na(anno_df$external_gene_name) | anno_df$external_gene_name == "",
  anno_df$ensembl_gene_id,
  anno_df$external_gene_name
)

logCPM_sym <- logCPM
rownames(logCPM_sym) <- gene_symbols

logCPM_sym_df <- as.data.frame(logCPM_sym)
logCPM_sym_df$Gene <- rownames(logCPM_sym_df)

logCPM_sym_collapsed <- logCPM_sym_df %>%
  group_by(Gene) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

logCPM_sym_collapsed <- as.data.frame(logCPM_sym_collapsed)
rownames(logCPM_sym_collapsed) <- logCPM_sym_collapsed$Gene
logCPM_sym_collapsed$Gene <- NULL

cat("\nExample rownames in logCPM_sym_collapsed:\n")
print(head(rownames(logCPM_sym_collapsed), 20))

# -------------------------------
# STEP 5 — IEG + Nur77 heatmap
# -------------------------------

ieg_genes <- c(
  "Zfp36", "Zfp36l1",
  "Pdcd1", "Irf4", "Il2ra", "Il2",
  "Icos", "Cd40lg", "Tfrc",
  "Myc", "Hif1a", "Il17a", "Ifng", "Il4",
  "Fos"  # c-Fos
)

nur77_gene <- "Nr4a1"

all_targets <- unique(c(ieg_genes, nur77_gene))

present_genes <- intersect(all_targets, rownames(logCPM_sym_collapsed))
missing_genes <- setdiff(all_targets, rownames(logCPM_sym_collapsed))

cat("\nNumber of IEG+Nur77 targets present:\n")
cat("  Present:", length(present_genes), "\n")
cat("  Missing:", length(missing_genes), "\n")

ieg_nur_present <- present_genes
mat_ieg_nur <- logCPM_sym_collapsed[ieg_nur_present, , drop = FALSE]

# Pretty group labels for heatmap
sample_info_mln$Group_pretty <- dplyr::recode(
  as.character(sample_info_mln$Group),
  "WT_PBS" = "WT +PBS",
  "WT_aGC" = "WT+αGalCer",
  "KO_PBS" = "CD1d-KO + PBS"
)

sample_info_mln$Group_pretty <- factor(
  sample_info_mln$Group_pretty,
  levels = c("WT +PBS", "WT+αGalCer", "CD1d-KO + PBS")
)

annotation_col <- data.frame(
  Group = sample_info_mln$Group_pretty
)
rownames(annotation_col) <- rownames(sample_info_mln)

ord <- order(annotation_col$Group)
annotation_col <- annotation_col[ord, , drop = FALSE]
mat_ieg_nur <- mat_ieg_nur[, rownames(annotation_col), drop = FALSE]

annotation_colors <- list(
  Group = c(
    "WT +PBS"        = "#a6cee3",
    "WT+αGalCer"     = "#fdbf6f",
    "CD1d-KO + PBS"  = "#cab2d6"
  )
)

pdf(file.path(out_dir, "MLN_IEG_Nur77_heatmap_prettyGroups.pdf"), width = 6, height = 7)
pheatmap(
  mat_ieg_nur,
  scale = "row",
  clustering_method = "complete",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col    = annotation_col,
  annotation_colors = annotation_colors,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "TCR IEGs + Nur77 in mLN",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
)
dev.off()

cat("\nDone. Outputs written to:\n", normalizePath(out_dir), "\n")
