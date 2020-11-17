####################################################################
# CCM bioinfo tutorial -- NGS analysis in R
# 2020.11.17
# Part 1: Biostrings in R
# Author: Huayun Hou (huayun.hou@sickkids.ca)
# Resource: 
#   https://bioconductor.org/packages/release/bioc/html/Biostrings.html
#   https://web.stanford.edu/class/bios221/labs/biostrings/lab_1_biostrings.html
#   https://girke.bioinformatics.ucr.edu/GEN242/mydoc_Rsequences_04.html
####################################################################

# package installation ----------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("Biostrings", "GenomicRanges", "GenomicFeatures", "rtracklayer", "org.Hs.eg.db","AnnotationHub","TxDb.Hsapiens.UCSC.hg19.knownGene","Organism.dplyr", "ShortRead","Rsubread","Gviz","mosaics"))

# set working directory
setwd("~/Documents/NGS_in_R_tutorial/") # set it to the downloaded "NGS_in_R_tutorial" directory




# Working with strings in base R -------------------------------

# creating some DNA sequences consists of TCAG, length between 10-20
set.seed(1234)
seqs <- sapply(1:3, function(x) paste(sample(c("T","C","A","G"), size = sample(20:30, size = 1), replace = TRUE), collapse = ""))

# string manipulations in base R
nchar(seqs) # length of the string
substr(seqs, start = 3, stop = 5) # subset string
grep(pattern = "CAC", seqs) # which sequence contains a pattern
gregexpr("AA", seqs) # find all matching patterns in all strings 
gsub("AA", "NN", seqs) # string substitution

?grep

# introduction to Biostrings --------------------------------------------------------------
# Biologically meaningful string manipulation support
library(Biostrings)

# The Xstring class
?XString

# XstringSet object holds multiple strings
?XStringSet

dnaseq <- DNAStringSet(seqs)

# Similar string manipulation methods as in base R
width(dnaseq) # length of the strings
subseq(dnaseq, start = 3, end = 5) # subset string
vmatchPattern("AA", dnaseq) # pattern matching
pmatch_result <- vmatchPattern("CAA", dnaseq, max.mismatch = 1) # more flexible options!

# Various methods specific to biological strings!
methods(class = "XStringSet")

reverseComplement(dnaseq) # get the reverseComplement sequences 
letterFrequency(dnaseq, "GC", as.prob = TRUE) # get letter frequencies
translate(dnaseq) # translate DNA to AA

# trim left/right flanking patterns  
dnaseq_rpattern <- DNAStringSet(paste0(dnaseq, "AATTA"))
dnaseq_rpattern
trimLRPatterns(Rpattern = "AATTA", subject = dnaseq_rpattern)

# create an RNAstring
rnaseq <- RNAStringSet(dnaseq)
rnaseq
alphabet(rnaseq)
alphabet(dnaseq)



