# NGS_in_R_tutorial
This repository contains the scripts and data used for the CCM bioinformatics tutorial series -- NGS data analysis in R

## Package installation

```R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("Biostrings", "GenomicRanges", "GenomicFeatures", "rtracklayer", "org.Hs.eg.db","AnnotationHub","TxDb.Hsapiens.UCSC.hg19.knownGene","EnsDb.Hsapiens.v75", "ShortRead","Rsubread","Gviz","mosaics"))
```
