#' Find chimeric reads
#' 
#' Assumes the the GAlignments object does not contain multimapping reads.
#' That is, read names that appear more than ones in the file are considered 
#' chimeras.
#' 
#'
#'@param bam A GAlignments object, must include names
#'@author Helen Lindsay
#'@export
findChimeras <- function(bam){
  chimera_idxs <- which((duplicated(names(bam)) | 
                         duplicated(names(bam), fromLast = TRUE)))
  chimera_idxs <- chimera_idxs[order(as.factor(names(bam)[chimera_idxs]))]
  return(chimera_idxs)
}