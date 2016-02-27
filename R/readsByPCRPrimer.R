#'@title Finds overlaps between aligned reads and PCR primers
#'@description Short reads amplified with PCR primers should start and end
#'at defined positions.  However, the ends of an aligned read may be clipped
#'as sequencing technologies are prone to making errors at the start and end.
#'readsByPCRPrimer extrapolates the genomic location of entire reads
#'from their aligned sections by adding clipped sections, then finds near
#'exact matches to a set of PCR primers.  Note that this is not always a good
#'assumption, and is misleading in the case of chimeric reads where sections
#'clipped in one part of a chimera are aligned in another.
#'@param bam A set of aligned reads
#'@param primers A set of ranges that the unclipped reads may map to
#'@param ... Additional arguments
#'@author Helen Lindsay
#'@rdname readsByPCRPrimer
#'@export
setGeneric("readsByPCRPrimer", function(bam, primers, ...) {
           standardGeneric("readsByPCRPrimer")})

#'@param tolerance Number of bases by which reads and primers may differ
#'at each end (Default: 0)
#'@param verbose Print number of full and partial matches (Default: TRUE)
#'@param ignore.strand Passed to \code{\link[GenomicAlignments]{findOverlaps}}
#'and \code{\link[GenomicRanges]{disjoin}}. Should strand be ignored when
#'finding overlaps.  (Default: TRUE)
#'@param allow.partial Should reads that do not match the PCR boundaries, but
#'map to a region covered by only one primer be considered matches?  (Default: TRUE)
#'@param chimera.idxs Indices of chimeric reads within the bam.  If specified,
#'chimeras overlapping multiple pcr primers will be removed.
#'@return A \code{\link[S4Vectors]{Hits}} object where "query" is the index with
#'respect to bam and "subject" is the index with respect to the primers.
#'@seealso \code{\link[GenomicRanges]{GRanges}}, \code{\link[GenomicAlignments]{GAlignments}}
#'@rdname readsByPCRPrimer
setMethod("readsByPCRPrimer", signature("GAlignments", "GRanges"),
          function(bam, primers, ..., tolerance=0, verbose = TRUE,
                   ignore.strand = TRUE, allow.partial = TRUE,
                   chimera.idxs = NULL){

            bamgr <- addClipped(bam)
            readsByPCRPrimer(bamgr, primers, tolerance = tolerance,
                   verbose = verbose, ignore.strand = ignore.strand,
                   allow.partial = allow.partial, chimera.idxs = chimera.idxs, ...)
        })


#'@rdname readsByPCRPrimer
setMethod("readsByPCRPrimer", signature("GRanges", "GRanges"),
          function(bam, primers, ..., tolerance=0, verbose = TRUE,
                   ignore.strand = TRUE, allow.partial = TRUE,
                   chimera.idxs = NULL){

            if (verbose){
              message("Finding overlaps between reads and primers\n")
            }
            hits_pcr <- findOverlaps(bam, primers, ignore.strand = ignore.strand,
                                     type = "equal", maxgap = tolerance)
            if (verbose){
              hits_pcr_l <- length(unique(queryHits(hits_pcr)))
              message(sprintf("%s from %s (%.2f%%) reads overlap a pcr region
                              almost exactly\n",
                          hits_pcr_l, length(bam), hits_pcr_l/length(bam)*100))
            }
            if (any(duplicated(queryHits(hits_pcr)))){
              warning("Cannot distinguish all pcr targets with this tolerance")
            }
            if (allow.partial){
              # find the unique regions of the pcr primers
              get_dups <- function(x) duplicated(x) | duplicated(x, fromLast = TRUE)

              disjoint <- disjoin(primers, ignore.strand = ignore.strand)
              dj_to_primer <- findOverlaps(disjoint, primers)
              is_dup <- get_dups(queryHits(dj_to_primer))
              dj_to_primer <- dj_to_primer[!is_dup]

              # find overlaps between the remaining reads and the unique regions
              remaining <- setdiff(c(1:length(bam)), queryHits(hits_pcr))
              rhits <- findOverlaps(bam[remaining], disjoint[queryHits(dj_to_primer)])
              is_dup <- get_dups(queryHits(rhits))

              if (verbose & any(is_dup)) {
                message(sprintf("Excluded %s reads which map to unique regions of
                                multiple primers\n",
                    length(unique(queryHits(rhits[is_dup])))))
              }
              rhits <- rhits[!is_dup]
              rhits_to_primer <- S4Vectors::remapHits(rhits,
                  Lnodes.remapping = factor(remaining,
                                            levels = c(1:length(bam))),
                  Rnodes.remapping = factor(subjectHits(dj_to_primer),
                                            levels = c(1:length(primers))))
              if (verbose){
                rhitsl <- length(unique(queryHits(rhits)))
                message(sprintf("Of the %s reads that do not exactly match a pcr region,
%s (%.2f%%) partially overlap exactly one pcr region\n\n",
                    (length(bam)-hits_pcr_l),
                    rhitsl, rhitsl/(length(bam)-hits_pcr_l)*100))
              }
              hits_pcr <- union(hits_pcr, rhits_to_primer)
            }
            if (!is.null(chimera.idxs)){
              hits_pcr <- rmMultiPCRChimera(names(bam), hits_pcr, chimera.idxs,
                                            verbose = verbose)
            }
            if (verbose){
              message(sprintf("Found matches for %s reads from %s\n\n",
                          length(hits_pcr), length(bam)))
            }
            hits_pcr
          })


#'@title Extrapolates mapping location from clipped, aligned reads
#'@description Extrapolates the mapping location of a read by assuming that
#'the clipped regions should map adjacent to the mapped locations.  This
#'is not always a good assumption, particularly in the case of chimeric reads!
#'@author Helen Lindsay
#'@rdname addClipped
setGeneric("addClipped", function(bam, ...) {
  standardGeneric("addClipped")})

#'@param bam A GAlignments object
#'@param ... additional arguments
#'@return A \code{\link[GenomicRanges]{GRanges}} representation of the extended
#'mapping locations
#'@rdname addClipped
setMethod("addClipped", signature("GAlignments"),
          function(bam, ...){
            l_clip <- as.numeric(gsub(".*[A-Z].*", 0, gsub("[HS].*","", cigar(bam))))
            r_clip <- as.numeric(gsub("^$", 0, gsub('.*M|[HS]$', "", cigar(bam))))
            bamgr <- GenomicRanges::GRanges(seqnames(bam),
                              IRanges(start(bam) - l_clip, end(bam) + r_clip))
            names(bamgr) <- names(bam)
            bamgr
          })



#'@title Remove chimeric reads overlapping multiple primers
#'@description Finds and removes sets of chimeric read alignments
#' that overlap more than one guide, i.e. that cannot be unambiguously
#' assigned to a single guide. 
#'@param readnames A set of read names, used for identifying chimeric read sets
#'@param pcrhits A mapping between indices of reads and a set of pcr primers
#'@param chimera_idxs location of chimeric reads within the bam
#'@param ... Additional arguments
#'@return pcrhits, with chimeric reads mapping to different primers omitted.
#'@author Helen Lindsay
#'@rdname rmMultiPCRChimera
setGeneric("rmMultiPCRChimera", function(readnames, pcrhits, chimera_idxs, ...) {
  standardGeneric("rmMultiPCRChimera")})


#'@param verbose Display information about the chimeras (Default: TRUE)
#'@rdname rmMultiPCRChimera
setMethod("rmMultiPCRChimera", signature("character", "Hits", "integer"),
          function(readnames, pcrhits, chimera_idxs, ..., verbose = TRUE){
            chs <- chimera_idxs[chimera_idxs %in% queryHits(pcrhits)]
            get_dups <- function(x) duplicated(x) | duplicated(x, fromLast = TRUE)
            nms <- readnames[chs]
            is_dup <- get_dups(nms)
            dnms <- nms[is_dup]
            ch_to_pcr <- match(chs[is_dup], queryHits(pcrhits))
            pcr <- subjectHits(pcrhits)[ch_to_pcr]
            pcr <- rle(paste(pcr, dnms, sep = "."))
            dnms <- rle(dnms)
            one_primer <- rep(dnms$lengths, dnms$lengths) == rep(pcr$lengths, pcr$lengths)
            if (verbose){
              message(sprintf("%s (%.2f%%) of the chimeric reads overlap a pcr primer.
  For %s (%2.f%% reads) of these, at least 2 segments of the chimera
overlap a primer.
    Of these, %s (%.2f%%) overlap different primers\n",
                   length(chs), length(chs)/length(chimera_idxs) * 100,
                   sum(is_dup), sum(is_dup)/length(chs)*100,
                   sum(!one_primer), length(unique(dnms[!one_primer]))))
            }
            pcrhits[-ch_to_pcr[!one_primer]]
          })
