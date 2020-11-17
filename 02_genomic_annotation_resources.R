####################################################################
# CCM bioinfo tutorial -- NGS analysis in R
# 2020.11.17
# Part 2: Genomic annotation resources in R
# Author: Huayun Hou (huayun.hou@sickkids.ca)
# Resource: 
#   https://www.bioconductor.org/packages/release/workflows/vignettes/annotation/inst/doc/Annotation_Resources.html
####################################################################


# AnnotationDbi: orgDb and TxDb -------------------------------------------

# OrgDb: Organism level, gene centric, stores mapping between identifiers and other information.
library(org.Hs.eg.db)

# Methods for AnnotationDb object
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

# TxDb: transcript oriented, maps between genes and transcripts, and their genomic intervals
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Exercise: use columns() and keytypes() to explore the txdb object
#           use select to obtain the transcript name and genomic locations of a few genes


# use select() txdb
some_genes <- keys(txdb, keytype = "GENEID")[100:105]
select(txdb, 
       keys = some_genes,
       keytype = "GENEID",
       columns = c("GENEID","TXNAME","TXCHROM","TXSTART","TXEND","TXSTRAND"))

# txdb specific methods using the GenomicFeatures package

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

# Given some genes, find their promoter regions 
some_genes <- keys(org.Hs.eg.db, keytype = "SYMBOL")[1:10]

gene_entrez <- select(org.Hs.eg.db, keys = some_genes, keytype = "SYMBOL", column = "ENTREZID")

gene_prom <- promoters(txdb, 
          upstream = 3000, 
          downstream = 3000, 
          filter = list(gene_id = gene_entrez$ENTREZID),
          columns = c("tx_name","gene_id"))



# match between gene entrez id and symbol
temp <- gene_entrez$SYMBOL
names(temp) <- gene_entrez$ENTREZID

# add to Grange object 
gene_prom$gene_symbol <- temp[unlist(gene_prom$gene_id)]

# Making TxDb object from gtf/gff files and more about TxDb objects, see:
# https://bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf

# In addition, EnsDb packages provides ensembl centric genes annotations (gene name, id, transcripts, location, biotype etc.) and can be access using the same select() method. 
# https://bioconductor.org/packages/release/bioc/html/ensembldb.html

# organism.dplyr ----------------------------------------------------------
# Since OrgDb and TxDb contain different information, it is bothersome to match between the two packages. Organism.dplyr makes a sqlite database of a organism from its corresponding orgDb and TxDb, queries both databases and merges the results, making it easier to get all information from "one resource". 

library(Organism.dplyr)

# see supported organisms that have both OrgDb and TxDb:
supportedOrganisms()

# initiate a src_organism object
src <- src_organism("TxDb.Hsapiens.UCSC.hg19.knownGene")

columns(src)

# Methods 

# 1.similar as in OrgDb and TxDb, select() can be used to retrieve information from an src_organism object

select(src, 
       keys = "TP53",
       keytype = "symbol",
       columns = c("alias","tx_name", "tx_start"))

# 2. methods from TxDb packages can be used with additional filters
supportedFilters(src)

# create filters
gr_filter <- GRangesFilter(GRanges("chr1:30000000-40000000"))
symbol_filter <- SymbolFilter("A", "startsWith")

promoters(src, 
          filter = ~gr_filter & symbol_filter)

transcripts(src, 
            filter = ~symbol %startsWith% "A" & gr_filter)

# Example:get promoter sequences for all transcripts on chr17

# tx_chrom filter does not work on "promoters()" directly, thus we get transcripts first, and use the promoters() function to get promoters. 
chr17_promoters <- promoters(transcripts(src,
                             filter = ~tx_chrom == "chr17",
                             columns = c("symbol")),
                            upstream = 3000,
                            downstream = 3000)

# 3. dplyr methods 
# count the number of transcripts for each gene on chr22
# see vignette for more information

tbl(src, "ranges_tx") %>% 
        filter(tx_chrom == "chr17" & !is.na(entrez)) %>% 
        left_join(tbl(src, "id"), by = "entrez") %>% 
        dplyr::select(symbol,ensembl, tx_name) %>% 
        distinct() %>% 
        group_by(symbol, ensembl) %>% 
        summarise(nTranscripts = n()) %>% 
        arrange(desc(nTranscripts))



# BSgenomes ---------------------------------------------------------------
# BSgenomes packages contain sequence data for a specific genome build of an organism
library(BSgenome)
available.genomes()
installed.genomes()

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library("BSgenome.Hsapiens.UCSC.hg19")

Hsapiens

# see all chromosomes
seqnames(Hsapiens)

# access sequence of one chromosome
chr17_seq <- Hsapiens[["chr17"]]

# convert to a XStringSet to save as a fasta file (we can only write XstringSet, not Xstring)
chr17_seq <- DNAStringSet(chr17_seq)
names(chr17_seq) <- "chr17"

# write chr17 sequence into a gzipped fasta file
writeXStringSet(chr17_seq, "hg19_chr17.fasta.gz", format = "fasta", compress = T)


# we can use Biostring functions on the sequence
# calculate GC content
letterFrequency(chr17_seq, "GC", as.prob = T)


# get promoter sequences
getSeq(Hsapiens, chr17_promoters)

# AnnotationHub -----------------------------------------------------------
# A client interface to abundant resources stored at the AnnotationHub web service 

library(AnnotationHub)
# initiate an AnnotationHub instance
ah <- AnnotationHub()

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






