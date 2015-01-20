#'@title Reverses the order of operations in a cigar string
#'@description For example, the string "20M5D15M" would become "15M5D20M"
#'@param cigar the cigar string.  
reverseCigar <- function(cigar){
  cigar.widths <- rev(strsplit(cigar, '[A-Z]')[[1]])
  cigar.ops <- rev(explodeCigarOps(cigar)[[1]])
  paste0(cigar.widths,cigar.ops, collapse = "")
}

#ALTERNATIVE IMPLEMENTATION TO TEST (VECTORISED)
#wdths <- explodeCigarOpLengths(cigs)
#ops <- explodeCigarOps(cigs)
#temp <- rev(relist(rev(paste0(unlist(wdths), unlist(ops))), rev(wdths)))
#result <- Map(paste, temp, collapse = "")
