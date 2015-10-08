#'@title Find chimeric reads
#'@description Find chimeric reads, assuming that the GAlignments object
#' does not contain multimapping reads. That is, read names that appear
#' more than ones in the file are considered chimeras.  Chimeric reads
#' are reads that cannot be mapped as a single, linear alignment.  Reads
#' from structual rearrangements such as inversions can be mapped as chimeras.
#' Note that the indices of all chimeric reads are returned, these are not
#' separated into individual chimeric sets.
#'@param bam  A GAlignments object, must include names
#'@param by.flag Can the chimeras be detected just using the supplementary
#'alignment flag?  (Default: FALSE).  If TRUE, detects supplementary alignments
#'and returns reads with the same name as a supplementary alignment (quicker).
#'If FALSE, all alignments with duplicated names are returned.
#'@author Helen Lindsay
#'@return A vector of indices of chimeric sequences within the original bam
#'@seealso \code{\link{plotChimeras}} for plotting chimeric alignment sets.
#'@export
#'@examples
#'bam_fname <- system.file("extdata", "gol_F1_clutch_2_embryo_4_s.bam",
#'                          package = "CrispRVariants")
#'bam <- GenomicAlignments::readGAlignments(bam_fname, use.names = TRUE)
#'chimera_indices <- findChimeras(bam)
#'chimeras <- bam[chimera_indices]
findChimeras <- function(bam, by.flag = FALSE){
  if (isTRUE(by.flag) & ! "flag" %in% names(mcols(bam))){
    stop("If 'by.cigar' is TRUE, bam must have a metadata column 'flag'")
  }
  if (isTRUE(by.flag)){
    suppl <- bitwAnd(mcols(bam)$flag, 2048)
    ch_names <- names(bam[suppl != 0])
    chimera_idxs <- which(names(bam) %in% ch_names)
  } else {
    chimera_idxs <- which((duplicated(names(bam)) |
                           duplicated(names(bam), fromLast = TRUE)))
  }
  chimera_idxs <- chimera_idxs[order(as.factor(names(bam)[chimera_idxs]))]
  chimera_idxs
}

