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

if(!dir.exists("index")){
  dir.create("index")  
}

# building index. For this tutorial, we are only building an index of hg19 chromosome 17, using the sequence we retrieved in the previous session
buildindex(basename="index/hg19_chr17_index",
           reference="data/hg19_chr17.fasta.gz", 
           gappedIndex = TRUE) # use gappedIndex to reduce computation resources needed


# Step3: alignment  -------------------------------------------------------

# create a directory for alignment outputs
if(!dir.exists("aln")){
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

# After alignment, we usually want to sort the bam files, create index for it so other programs can efficiently access the file content, and sometimes we might want to further filter the file, for example, keep alignments that are above a certain mapping quality, or keeping uniquely mapped reads only, for downstream analysis.
# We can further manipulate bam files using Rsamtools
library(Rsamtools)

bam_files <- dir("aln", pattern = "*bam$", full.names = T)
sorted_bam_prefix <- gsub("bam", "sorted", bam_files)
sorted_bam_files <- paste0(sorted_bam_prefix, ".bam")
for (i in 1:3){
  sortBam(bam_files[[i]],
        destination = sorted_bam_prefix[[i]])
}

indexBam(sorted_bam_files)

library(GenomicAlignments)
bamfile <- "aln/RelA_chip_rep1.bam"
alignment <- readGAlignments(bamfile)

flags <- scanBamFlag(isSecondaryAlignment = F)
alignment <- readGAlignments(bamfile, param = ScanBamParam(flag = flags, what = "seq"))

# Step4: peak calling -----------------------------------------------------
library(mosaics)

# create a directory for alignment outputs
if(!dir.exists("peaks")){
  dir.create("peaks")  
}

# specify parameters
fragment_length <- 200
bin_size <- 200

# We will use one ChIP experiment file and one corresponding input (control) file
input_file <- sorted_bam_files[[1]]
chip_file <- sorted_bam_files[[2]]

# First, we construct bins across regions covered in our ChIP and input files
constructBins(infile = chip_file, 
              fileFormat = "bam",
              outfileLoc = "peaks",
              fragLen = fragment_length,
              binSize = bin_size)

constructBins(infile = input_file, 
              fileFormat = "bam",
              outfileLoc = "peaks",
              fragLen = fragment_length,
              binSize = bin_size)

# get the names of these bin files
chip_bin_count_file <- paste0("peaks/", basename(chip_file),"_fragL",fragment_length,"_bin",bin_size,".txt")
input_bin_count_file <- paste0("peaks/", basename(input_file),"_fragL",fragment_length,"_bin",bin_size,".txt") 

# count numbers of reads mapping to each bin
bin_count <- readBins(type = c("chip", "input"), fileName = c(chip_bin_count_file, input_bin_count_file))

# fit the "mosaics" model
fit <- mosaicsFit(bin_count, analysisType = "IO")

# finding peaks (with default parameters)
peaks <- mosaicsPeak(fit)

# access peaks using the print() method
head(print(peaks))
peaks_df <- print(peaks)

# convert peaks from data frame to GRanges
# first we need to rename the chromosome, start, and end columns
names(peaks_df) <- c("seqnames", "start", "end")
# convert data frame to GRanges
peaks_gr <- makeGRangesFromDataFrame(peaks_df)

# export peaks as bed files
export(peaks, type = "bed", filename = paste0("peaks/", basename(chip_file), "_peaks.bed"))



# Step5: visualization -----------------------------------------------------------
# Our aim is to visualize the genomic location around a gene, showing the gene model, ChIP-seq signal (read coverage), and peaks called
library(Gviz)
# first we need to get the genomic coordinate we want to visualize. We are focusing on  gene CCL2, and 3kb upstream of its promoter.
gene_range <- transcripts(src, 
                          filter = ~symbol == "CCL2")
upstream_3k <- flank(gene_range, width = 3000)
plot_range <- range(c(gene_range, upstream_3k))
seqlevels(plot_range) <- "chr17"

# buiding tracks
# build an annotation track to show peak coordinates 
atrack <- AnnotationTrack(subsetByOverlaps(peaks_gr, plot_range), 
                          name = "Peaks",
                          rotation.title = 3)
# plot the track using "plotTracks" function
plotTracks(atrack)

# generate a genomic axis
gtrack <- GenomeAxisTrack()

# generate a gene model track from TxDb 
txTr <- GeneRegionTrack(txdb, chromosome = as.character(seqnames(plot_range)), 
                        start = start(plot_range),  end = end(plot_range), 
                        name = "gene model", 
                        transcriptAnnotation = "symbol",
                        rotation = 90)

# generate alignment track for the ChIP sample
alnTr <- AlignmentsTrack(sorted_bam_files[[2]], chromosome = as.character(seqnames(plot_range)), 
                         start = start(plot_range),  end = end(plot_range), fill = "orange", name = "RelA_ChIP")
# generate alignment track for the Input sample
alnTr2 <- AlignmentsTrack(sorted_bam_files[[1]], chromosome = as.character(seqnames(plot_range)), 
                         start = start(plot_range),  end = end(plot_range), name = "Input")

# plot all tracks together 
plotTracks(list(gtrack, atrack, txTr, alnTr, alnTr2), type = "coverage")



# An alternative package "ggbio" can be used to plot gene symbols, bam file coverage, and variants. The resulting plots can be easily combined with other ggplot objects for customized plot generation.

# Some example code
# p_gene <- autoplot(txdb, which = plot_range)
# 
# p_signal <- autoplot(sorted_bam_files[[1]], which = plot_range)
# 
# tks <- tracks(chip_signal = p_signal,
#               gene_model = p_gene)

