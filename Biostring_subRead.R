# package installation ----------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("systemPipeR")
BiocManager::install("systemPipeRdata")
source("https://bioconductor.org/biocLite.R")
BiocManager::install(c("Biostrings", "GenomicRanges", "GenomicRanges", "rtracklayer", "systemPipeR", "seqLogo", "ShortRead"))
BiocManager::install("Rsubread")


# load libraries
library(ShortRead)
library(Rsubread)

# Working with strings in base R or with Biostrings -------------------------------

# creating some DNA sequences consists of TCAG, length between 10-20
seqs <- sapply(1:3, function(x) paste(sample(c("T","C","A","G"), size = sample(10:20, size = 1), replace = TRUE), collapse = ""))

# string manipulations in base R
nchar(seqs) # length of the string
substr(seqs, start = 3, stop = 5) # subset string
grep(pattern = "CAC", seqs) # which sequence contains a pattern
gregexpr("AA", seqs) # find all matching patterns in all strings 
gsub("AA", "NN", seqs) # string substitution

# Biologically meaningful string manipulation support
library(Biostrings)
dnaseq <- DNAStringSet(seqs)

# set working directory
setwd("~/Documents/tutorial/NGS_in_R/")

fastq_files <- dir("data", "*fastq.gz", full.names = T)

# count number of reads in a fastq file 
countFastq(c(fastq_file1, fastq_file2))

# If the fastq file is small, we can read it into memory 
fq <- readFastq(fastq_file)

# access reads
sread(fq)[10:20]

qa_summary <- qa(fastq_files, type = "fastq")
browseURL(report(qa_summary))


# building index 
buildindex(basename="mm10_chr19_index",reference="data/chr19.fa")

# create a directory for outputs
if(!file.exists("aln")){
  dir.create("aln")  
}


# perform alignment 
alignment_log <- align(index = "mm10_chr19_index",
      readfile1 = fastq_files,
      type = "dna",
      output_file = paste0("aln/",gsub("fastq.gz", "bam", basename(fastq_files))), 
      unique = FALSE,
      nBestLocations = 5)


