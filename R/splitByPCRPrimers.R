#'@title Finds overlaps between aligned reads and PCR primers
#'@description Short reads amplified with PCR primers should start and end
#'at defined positions.  However, the ends of an aligned read may be clipped
#'as sequencing technologies are prone to making errors at the start and end.
#'splitByPCRPrimers extrapolates the genomic location of entire reads 
#'from their aligned sections by adding clipped sections, then finds near 
#'exact matches to a set of PCR primers.
#'@param bam A set of aligned reads
#'@param primer_ranges A set of ranges that the unclipped reads may map to 
#'@author Helen Lindsay
#'@rdname splitByPCRPrimers
#'@export
setGeneric("splitByPCRPrimers", function(bam, primers,...) {
           standardGeneric("splitByPCRPrimers")})

#'@param bam A GAlignments object 
#'@param primer_ranges A GRanges object 
#'@param tolerance Number of bases by which reads and primers may differ 
#'at each end (Default: 0)
#'@param verbose Print number of full and partial matches (Default: TRUE)
#'@return A \code{\link[IRanges]{Hits}} object where "query" is the index with 
#'respect to bam and "subject" is the index with respect to the primers.
#'@seealso \code{\link[GenomicRanges]{GRanges}}, \code{\link[GenomicAlignments]{GAlignments}}
#'@rdname splitByPCRPrimers
setMethod("splitByPCRPrimers", signature("GAlignments", "GRanges"),
          function(bam, primers, ..., tolerance=0, verbose = TRUE){
           l_clip <- as.numeric(gsub(".*[A-Z].*", 0, gsub("[HS].*","", cigar(bam))))
           r_clip <- as.numeric(gsub("^$", 0, gsub('.*M|[HS]$', "", cigar(bam)))) 
            bamgr <- GRanges(seqnames(bam), IRanges(start(bam) - l_clip, end(bam) + r_clip))
            
            hits_pcr <- findOverlaps(bamgr, primers, ignore.strand = TRUE, 
                                     type = "equal", maxgap = tolerance)
            
            hits_pcr_l <- length(unique(hits_pcr@queryHits))
            if (verbose == TRUE){
              cat(sprintf("%s from %s (%.2f%%) reads overlap a pcr region\n", 
                          hits_pcr_l, length(bam), hits_pcr_l/length(bam)*100))
              remaining <- setdiff(c(1:length(bam)), hits_pcr@queryHits) 
              rhits <- findOverlaps(bamgr[remaining], pcr_targets)  
              rhitsl <- length(unique(rhits@queryHits))   
              cat(sprintf("Of the remaining %s reads, %s (%.2f%%) partially overlap a pcr region\n\n", 
                          (length(bam)-hits_pcr_l),rhitsl, rhitsl/(length(bam)-hits_pcr_l)*100))
            }
            if (any(duplicated(hits_pcr@queryHits))){
              warning("Cannot distinguish pcr targets with this tolerance")
            } 
           return(hits_pcr)
        })