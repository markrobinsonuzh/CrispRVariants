#'@title Removes reads from a bam file 
#'@description Returns a GAlignments excluding reads based on either name and/or location
#'@param bam a GAlignments object
#'@param exclude.ranges Regions to exclude, as \code{\link[GenomicRanges]{GRanges}}.
#'@param exclude.names A character vector of alignments names to exclude
#'@rdname excludeFromBam
#'@return The bam minus the excluded regions
#'@author Helen Lindsay
excludeFromBam <- function(bam, exclude.ranges = GRanges(), exclude.names = NA){
  if (length(exclude.ranges) > 0) bam <- excludeFromBamByRange(bam, exclude.ranges)
  if (! is.na(exclude.names)) bam <- excludeFromBamByName(bam, exclude.names)
  return(bam)
}

excludeFromBamByName <- function(bam, exclude.names){
  excluden <- which(names(bam) %in% exclude.names)     
  bam <- bam[setdiff(seq_along(bam), excluden)]
}

excludeFromBamByRange <- function(bam, exclude.ranges){
  fo <- queryHits(findOverlaps(bam, exclude.ranges))
  return(bam[setdiff(seq_along(bam), fo)])
}