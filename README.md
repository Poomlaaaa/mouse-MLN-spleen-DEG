# Study Overview
Gene expression comparison between Mesenteric Lymph Node (MLN) and Spleen iNKT cells.
The objective is to investigate transcriptional differences between iNKT cells isolated from the mesenteric lymph nodes (MLN) and the spleen of wild-type mice, with a focus on the expression of **Nur77 (Nr4a1)** and **endoplasmic reticulum (ER) stress-related genes**.
The goal is to determine whether MLN iNKT cells exhibit a higher baseline stress or activation state compared to their splenic counterparts.

## Experimental Design
- **Groups:** Naïve Spleen (control, treated as *healthy*) vs. Naïve MLN (treated as *unhealthy*).
- **Biological question:** Identify genes that are upregulated or downregulated between MLN and spleen iNKT cells, and perform pathway enrichment analysis based on these findings.
- **Data source:** Publicly available RNA-seq datasets from the NCBI Gene Expression Omnibus (GEO) and Sequence Read Archive (SRA).
- **Focus:** Stress-related, tissue-adapted, or regulatory gene expression differences under homeostatic (non-stimulated) conditions.

## Analysis Pipeline
1. **Input:** STAR `ReadsPerGene` counts (from GEO/SRA data).
2. **Differential expression analysis:** `edgeR`.
3. **Gene annotation:** via `biomaRt`.
4. **Filtering:** Focus on ER stress-related genes.
5. **Pathway analysis:** GO enrichment using `clusterProfiler`.
6. **Visualization:** Heatmap of selected gene expression patterns.
7. **Output:** CSV tables of DEGs, annotations, ER stress gene lists, and enrichment results.

## Required R Packages

```r
install.packages(c("dplyr", "pheatmap"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("edgeR", "clusterProfiler", "org.Mm.eg.db", "biomaRt"))
