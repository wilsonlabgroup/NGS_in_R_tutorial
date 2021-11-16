####################################################################
# CCM bioinfo tutorial -- NGS analysis in R
# 2020.11.17
# Part 3: Genomic range operations in R
# Author: Huayun Hou (huayun.hou@sickkids.ca)
# Resource: 
#   - https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html
####################################################################

# GRanges: object that represent genomic ranges
my_promoters <- chr17_promoters[1:5]
class(my_promoters)
my_promoters

# GRanges is built upon IRanges: a fundemental structure in bioconductor for defining integer ranges. IRanges contain start, end, and width
class(pmatch_result[[1]])
pmatch_result[[1]]

# While IRanges specify ranges in general, GRanges add chromosome and strand information to address genomic locations 
# 1-based: coordinate starts counting at 1
# left-most: regions on the minus strands are defined to start at the left-most coordinate
# closed: region includes both the start and end position; a region with the same start and end has a width of 1 


# methods on GRanges:
# Subset:
my_promoters[1:3] # subset by index, only one dimension 
my_promoters[2:5, "symbol"] # subset by metadata column name
names(my_promoters) <- my_promoters$tx_name # naming GRanges
my_promoters["ENST00000343572"] # subset by name 

# access ranges elements
seqnames(my_promoters)
start(my_promoters)
strand(my_promoters)


# get metadata (columns beyond seqnames, ranges, and strand)
mcols(my_promoters) # it is a data frame
# change metadata; keeping only "tx_name" and "symbol" columns
mcols(my_promoters) <- mcols(my_promoters)[, c("tx_name","symbol")]
# metadata columns can also be accessed using $
my_promoters$tx_name

# explore GRangse
width(my_promoters)
length(my_promoters) # numbers of ranges in an GRanges object
seqlevels(my_promoters) # show all chromosomes stored in the object



# Useful interval manipulation methods
# 1. Intra-range methods: methods that apply to every range in a GRanges object

?"intra-range-methods"
flank(my_promoters, width = 10) # get N bp upstream
shift(my_promoters, shift = -10) # shift the region by N bp
resize(my_promoters, width = 100, fix = "start") # resize 


# 2. inter-range methods: comparing between ranges within a single GRanges object
reduce(my_promoters) # collapse overlapping ranges 
gaps(my_promoters) # get the gaps between ranges 
coverage(my_promoters) # get the coverage/degress of overlap along the sequence

# 3. between-range methods: calculates the relationships between two GRanges objects, for example, overlaps

# use findOverlaps to find overlaps between a query set (first) and a subject set (second)

# change seqlevels to "UCSC style"
seqlevelsStyle(chr17_promoters) <- "UCSC"

# find overlaps
ol <- findOverlaps(peak_set1, chr17_promoters)
ol # a Hits object containing the index pairing for the query set and the subject set 

# get the original overlapping regions

# use from() to get query set index, use to() to get subject set index
# or transform Hits objects to a data frame and subset 

peaks_ol_promoters <- peak_set1[from(ol)] # to(ol)
peaks_ol_promoters <- peak_set1[as.data.frame(ol)$queryHits]

# get metadata from "prmoters" and add them to overlapping peaks
promoters_ol_peaks <- chr17_promoters[to(ol)]
mcols(peaks_ol_promoters) <- cbind(mcols(peaks_ol_promoters),
                                   mcols(promoters_ol_peaks))

# Other overlap methods
countOverlaps(peak_set1, chr17_promoters) # count the number of subject overlaps for each query range
subsetByOverlaps(peak_set1, chr17_promoters) # returns only ranges in the query set that overlaps the subject set


# nearest methods
nearest(peak_set1, chr17_promoters) # for each range in the query set, finds the nearest entry in subject set 
distanceToNearest(peak_set1, chr17_promoters) # returns a Hits object with index from each query range and its nearest subject range, as well as the distance between the two features 

# see more methods
?"nearest-methods"




