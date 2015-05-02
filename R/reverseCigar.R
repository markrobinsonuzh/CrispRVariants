#'@title Reverses the order of operations in a cigar string
#'@description For example, the string "20M5D15M" would become "15M5D20M"
#'@param cigar the cigar string.  
#'@return The reversed cigar string
reverseCigar <- function(cigar){
  cigar.widths <- rev(strsplit(cigar, '[A-Z]')[[1]])
  cigar.ops <- rev(explodeCigarOps(cigar)[[1]])
  paste0(cigar.widths,cigar.ops, collapse = "")
}

#'@title Vectorised version of reverseCigar
#'@description For example, the string "20M5D15M" would become "15M5D20M"
#'reverseCigarv will replace reverseCigar after testing
#'@param cigars the cigar string.  
#'@return The reversed cigar string
reverseCigarv <- function(cigars){
  wdths <- explodeCigarOpLengths(cigars)
  ops <- explodeCigarOps(cigars)
  temp <- rev(relist(rev(paste0(unlist(wdths), unlist(ops))), rev(wdths)))
  result <- as.character(Map(paste, temp, collapse = ""))
}
