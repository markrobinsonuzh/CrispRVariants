#'Removes reads from a bam file 
#'
#'Returns a GAlignments excluding reads based on either name and/or location
#'
#'@param bam a GAlignments object
#'@param exclude_ranges
#'@param exclude_names
#'@author Helen Lindsay
#'@export
excludeFromBam <- function(bam, exclude_ranges = GRanges(), exclude_names = NA){
  if (length(exclude_ranges) > 0) bam <- excludeFromBamByRange(bam, exclude_ranges)
  if (! is.na(exclude_names)) bam <- excludeFromBamByName(bam, exclude_names)
  return(bam)
}

#'Removes reads from a bam file 
#'
#'Returns a GAlignments excluding reads based on name
#'
#'@param bam A GAlignments object
#'@param exclude_names a character vector of sequence names to exclude
#'@author Helen Lindsay
excludeFromBamByName <- function(bam, exclude_names){
  excluden <- which(names(bam) %in% exclude_names)     
  bam <- bam[setdiff(seq_along(bam), excluden)]
}

#'Removes reads from a bam file 
#'
#'Returns a GAlignments excluding reads based on location
#'
#'@param bam a GAlignments object
#'@param exclude_ranges A GRanges object containing regions to exclude 
#'@author Helen Lindsay
excludeFromBamByRange <- function(bam, exclude_ranges){
  fo <- findOverlaps(bam, exclude_ranges)@queryHits
  return(bam[setdiff(seq_along(bam), fo)])
}