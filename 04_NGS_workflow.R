####################################################################
# CCM bioinfo tutorial -- NGS analysis in R
# 2020.11.17
# Part 4: A minimal example of an NGS workflow and visualization entirely in R
# Author: Huayun Hou (huayun.hou@sickkids.ca)
# Resource: 
#   
####################################################################


# load libraries ----------------------------------------------------------
library(ShortRead) # package for inspection and quality control of fastq files
library(Rsubread) # package to perform alignment 


# For this tutorial, we are using ChIP-seq data of RELA, which encodes a core subunit of the nuclear factor NF-kB.
# Two biological replicates were performed. Only a small subset of the reads are used here as an example.

# Examine fastq files -----------------------------------------------------

# list fastq files in the data folder
fastq_files <- dir("data", "*fastq.gz", full.names = T)

# Some fastq file-handling functions from package "ShortRead"

# count number of reads in a fastq file. Here, we are not reading the fastq files into memory
countFastq(fastq_files)

# If the fastq file is small enough, we can actually read it into memory 
fq <- readFastq(fastq_files[[1]])

class(fq)
?"ShortRead"

# access fastq read information 
ShortRead::id(fq)[100:101] # read IDs
sread(fq)[100:101] # read sequences
quality(fq)[100:101] # per base quality
as(quality(fq), "matrix")[100:101, 1:40] # convert quality encoding to numeric quality scores



#  Step1: quality assessment  ---------------------------------------------
# Most commonly used fastq quality control tool is fastqc: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# Here, we used "ShortReads" package for fastq quality assessment. There are other R packages such as "Rqc" to perform quality control as well.

qa_summary <- qa(fastq_files, type = "fastq")
browseURL(report(qa_summary))

# NOTE: adapter and quality trimming: as we do not see high adapter content in the reads and read quality is high, we are skipping the trimming step here. 
# see "filtering and trimming" from ShorRead vignette for more information: http://bioconductor.org/packages/release/bioc/vignettes/ShortRead/inst/doc/Overview.pdf 



# Step2: generating genome index ------------------------------------------

# Before alignment, we need to generate an index for the corresponding genome. Each aligner requires a different index. We create a directory to store the index files

if(!file.exists("index")){
  dir.create("index")  
}

# building index. For this tutorial, we are only building an index of hg19 chromosome 17, using the sequence we retrieved in the previous session
buildindex(basename="index/hg19_chr17_index",
           reference="hg19_chr17.fasta.gz", 
           gappedIndex = TRUE) # use gappedIndex to reduce computation resources needed


# Step3: alignment  -------------------------------------------------------

# create a directory for alignment outputs
if(!file.exists("aln")){
  dir.create("aln")  
}


# perform alignment using Rsubread
alignment_log <- align(index = "index/hg19_chr17_index",
                       readfile1 = fastq_files, # for paired read, "readfile1" will be read1 from sequncer, and needs to be paried with read2 files
                       type = "dna",
                       output_file = paste0("aln/",gsub("fastq.gz", "bam", basename(fastq_files))), 
                       unique = FALSE,
                       nBestLocations = 5)

# examine alignment summary
View(alignment_log)

