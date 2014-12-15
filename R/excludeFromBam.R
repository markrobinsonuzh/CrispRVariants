#'Removes reads from a bam file 
#'
#'Returns a GAlignments excluding reads based on either name and/or location
#'
#'@param bam a GAlignments object
#'@param exclude.ranges
#'@param exclude.names
#'@author Helen Lindsay
excludeFromBam <- function(bam, exclude.ranges = GRanges(), exclude.names = NA){
  if (length(exclude.ranges) > 0) bam <- excludeFromBamByRange(bam, exclude.ranges)
  if (! is.na(exclude.names)) bam <- excludeFromBamByName(bam, exclude.names)
  return(bam)
}

#'Removes reads from a bam file 
#'
#'Returns a GAlignments excluding reads based on name
#'
#'@param bam A GAlignments object
#'@param exclude.names a character vector of sequence names to exclude
#'@author Helen Lindsay
excludeFromBamByName <- function(bam, exclude.names){
  excluden <- which(names(bam) %in% exclude.names)     
  bam <- bam[setdiff(seq_along(bam), excluden)]
}

#'Removes reads from a bam file 
#'
#'Returns a GAlignments excluding reads based on location
#'
#'@param bam a GAlignments object
#'@param exclude.ranges A GRanges object containing regions to exclude 
#'@author Helen Lindsay
excludeFromBamByRange <- function(bam, exclude.ranges){
  fo <- findOverlaps(bam, exclude.ranges)@queryHits
  return(bam[setdiff(seq_along(bam), fo)])
}