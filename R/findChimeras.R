#'@title Find chimeric reads
#'@description Find chimeric reads, assuming that the GAlignments object
#' does not contain multimapping reads. That is, read names that appear
#' more than ones in the file are considered chimeras.  Chimeric reads
#' are reads that cannot be mapped as a single, linear alignment.  Reads 
#' from structual rearrangements such as inversions can be mapped as chimeras.
#' Note that the indices of all chimeric reads are returned, these are not 
#' separated into individual chimeric sets.
#'@param bam  A GAlignments object, must include names
#'@author Helen Lindsay
#'@return A vector of indices of chimeric sequences within the original bam
#'@seealso \code{\link{plotChimeras}} for plotting chimeric alignment sets.
#'@export
#'@examples
#'bam_fname <- system.file("extdata", "gol_F1_clutch_2_embryo_4_s.bam",
#'                          package = "crispRvariants")
#'bam <- GenomicAlignments::readGAlignments(bam_fname, use.names = TRUE)
#'chimera_indices <- findChimeras(bam)
#'chimeras <- bam[chimera_indices]
findChimeras <- function(bam){
  chimera_idxs <- which((duplicated(names(bam)) | 
                         duplicated(names(bam), fromLast = TRUE)))
  chimera_idxs <- chimera_idxs[order(as.factor(names(bam)[chimera_idxs]))]
  return(chimera_idxs)
}