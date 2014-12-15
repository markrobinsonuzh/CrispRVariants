#'@title Reverses the order of operations in a cigar string
#'@description For example, the string "20M5D15M" would become "15M5D20M"
#'@param cigar the cigar string.  
reverseCigar <- function(cigar){
  cigar.widths <- rev(strsplit(cigar, '[A-Z]')[[1]])
  cigar.ops <- rev(explodeCigarOps(cigar)[[1]])
  paste0(cigar.widths,cigar.ops, collapse = "")
}