####################################################################
# CCM bioinfo tutorial -- NGS analysis in R
# 2021.11.15
# Part 2: Genomic annotation resources in R
# Author: Huayun Hou (huayun.hou@sickkids.ca)
# Resource: 
#   https://www.bioconductor.org/packages/release/workflows/vignettes/annotation/inst/doc/Annotation_Resources.html
# BioC workshop: http://biocworkshops2019.bioconductor.org.s3-website-us-east-1.amazonaws.com/page/Bioc2019Anno__AnnotationWorkshop/
####################################################################


# AnnotationDbi: orgDb and TxDb -------------------------------------------

# core functions: columns(), keytypes(), keys(), select()

# OrgDb: Organism level, gene centric, stores mapping between identifiers and other information.
library(org.Hs.eg.db)

# Methods for AnnotationDbi object
columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
head(keys(org.Hs.eg.db)) # entrez id 
head(keys(org.Hs.eg.db, keytype = "SYMBOL")) # get keys by keytype

# retrieve information for a gene 
gene_info <- select(org.Hs.eg.db, 
                    keys = "TP53", 
                    keytype = "SYMBOL", 
                    columns = c("SYMBOL","ENSEMBL","GO"))

# use mapIds() to extract only one column of data. Returns a named vector. 
some_genes <- keys(org.Hs.eg.db, keytype = "ENSEMBL")[1:20]
ensembl_to_symbol <- mapIds(org.Hs.eg.db, 
       keys = some_genes,
       keytype = "ENSEMBL",
       column = "SYMBOL")

# when there are multiple values to be returned? Use multiVals to determine behavior
ensembl_to_symbol_multi <- mapIds(org.Hs.eg.db, 
                            keys = some_genes,
                            keytype = "ENSEMBL",
                            column = "GO",
                            multiVals = "list")

# TxDb and EnsDb: transcript oriented, maps between genes and transcripts, and their genomic intervals
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# use columns() and keytypes() to explore the txdb object
# use select to obtain the transcript name and genomic locations of a few genes

# use select() txdb
set.seed(1234)
some_genes <- keys(txdb, keytype = "GENEID")[sample(1:23000, 10)]
select(txdb, 
       keys = some_genes,
       keytype = "GENEID",
       columns = c("GENEID","TXNAME","TXCHROM","TXSTART","TXEND","TXSTRAND"))

# txdb specific methods using the GenomicFeatures package
# https://bioconductor.org/packages/devel/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf
# TxDb objects can be made from gtf/gff files 

transcripts(txdb)[1:10]
exonsBy(txdb, by = "gene")[1:2]
promoters(txdb, 
          upstream = 3000, 
          downstream = 3000, 
          use.names = T, 
          filter = list(tx_chrom = "chr19"))[1:10]

# other useful methods
fiveUTRsByTranscript(txdb)[1:10]
threeUTRsByTranscript(txdb)
intronsByTranscript(txdb)

# EnsDb is similar to Txdb, but with ensembl-oriented annotations and additional functions
# Learn more at: https://bioconductor.org/packages/release/bioc/html/ensembldb.html

library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

edb

# With AnnotationDbi methods
columns(edb)
keytypes(edb)

some_genes <- keys(edb, keytype = "SYMBOL")[sample(1:20000, 10)]

gene_info <- select(edb, 
                    keytype = "SYMBOL",
                    keys = some_genes, 
                    columns = c("SYMBOL","GENEID","TXBIOTYPE","SEQNAME", "GENESEQSTART","GENESEQEND"))

# EnsDb specific:
# Filters
supportedFilters(edb)

# Filter with a "filter" object
tx <- transcripts(edb, filter = GeneNameFilter("TP53"))

# Filter with a formula: ~ <field> <condition> <value>
tx <- transcripts(edb, filter = ~gene_name == "TP53")

tx_pc <- transcripts(edb, filter = (~gene_name == "TP53" & tx_biotype == "protein_coding"))

# Specify columns to return from extractor methods ("genes","transcripts", "exons")

# get all promoter regions on chr17
chr17_promoters <- promoters(transcripts(edb, 
                                         filter = ~seq_name == 17 & tx_biotype == "protein_coding",
                                         columns = c("seq_name","tx_seq_start", "tx_seq_end", "tx_name","symbol","protein_id")))


# Other useful functions include converting between genome, transcript, and protein coordinate 

# BSgenomes ---------------------------------------------------------------
# BSgenomes packages contain sequence data for a specific genome build of an organism
library(BSgenome)
available.genomes()
installed.genomes()

# install a BSgenome
BiocManager::install("BSgenome.Celegans.UCSC.ce11")
library("BSgenome.Celegans.UCSC.ce11")

Celegans

# see all chromosomes
seqnames(Celegans)

# access sequence of one chromosome, the result is a DNAString object
chrI_seq <- Celegans[["chrI"]]

# write chr17 sequence into a gzipped fasta file
# convert to a XStringSet to save as a fasta file (we can only write XstringSet, not Xstring)
chrI_seq <- DNAStringSet(chrI_seq)
names(chrI_seq) <- "chrI"

# output
writeXStringSet(chrI_seq, "ce11_chrI.fasta.gz", format = "fasta", compress = T)

# we can use Biostring functions on the sequence
# calculate GC content
letterFrequency(chrI_seq, "GC", as.prob = T)

# get sequences
example_region <- GRanges(seqnames = "chrI", IRanges(100,200), strand = "-")
getSeq(Celegans, example_region)

# AnnotationHub -----------------------------------------------------------
# A client interface to abundant resources stored at the AnnotationHub web service 

library(AnnotationHub)
# initiate an AnnotationHub instance
ah <- AnnotationHub()
removeCache(ah)
ah

# use [] to get information about a specific dataset
ah[1000]


# explore 
unique(ah$dataprovider)

# query the database for OrgDb objects
all_species <- query(ah, "OrgDb")
# select data for human
all_species[grepl("Homo sapiens", all_species$species)]
# retrieve human orgDb object
# human_orgdb <- all_species[["AH84122"]]                 
            
# obtain data from epigenomeRoadMap 
epi <- query(ah, "epigenome")
# what types of files are available
table(epi$sourcetype)
h3k27ac_brain_peak <- query(ah, c("epigenome", "H3K27ac", "brain", "broadPeak"))       

# Subsetting annotationHub object using 'subset()' to find target datasets
target_ids <- subset(h3k27ac_brain_peak, grepl("E067|E068", title))$ah_id

# use [[]] to retrive data 
peak_set1 <- h3k27ac_brain_peak[[target_ids[1]]]
peak_set2 <- h3k27ac_brain_peak[[target_ids[2]]]

# or use display
display(h3k27ac_brain_peak)






