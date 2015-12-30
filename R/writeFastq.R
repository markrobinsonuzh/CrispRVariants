#'@title Append a sequence to a fastq file
#'
#'Used by abifToFastq
#'
#'@param outf Name of fastq file to append sequence
#'@param vals A list containing entries named "seq" (sequence) and
#'"quals" (quality scores, in ASCII format)
#'@param allow_spaces Should spaces in the sequence name be
#'substituted with underscores?  TRUE or FALSE
#'@return None.  The sequences in "vals" are written to outf
#'@author Helen Lindsay
writeFastq <- function(outf, vals, allow_spaces = FALSE){
  nms <- ifelse(allow_spaces, vals$seqname, gsub(" ", "_", vals$seqname))
  o <- file(outf, "a")
  seqname <- sprintf("@%s", nms)
  qualname <- sprintf("+%s", nms)
  writeLines(c(seqname, vals$seq, qualname, vals$quals) , o)
  close(o)
  return()
}