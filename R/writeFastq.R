#'@title Append a sequence to a fastq file
#'@description  Used by abifToFastq to write sanger sequences
#'to fastq format  As abifToFastq appends output to files,
#' writeFastq checks that sequence names are unique.  This
#' function is faster with checking switched off.
#'@param outf Name of fastq file to append sequence
#'@param vals A list containing entries named "seq" (sequence) and
#'"quals" (quality scores, in ASCII format)
#'@param allow_spaces Should spaces in the sequence name be
#'substituted with underscores?  TRUE or FALSE
#'@param check Check whether reads with the same name already
#'exist in the output fastq.  (Default: TRUE)
#'@return None.  The sequences in "vals" are written to outf
#'@author Helen Lindsay
writeFastq <- function(outf, vals, allow_spaces = FALSE, check = TRUE){
  nms <- ifelse(allow_spaces, vals$seqname, gsub(" ", "_", vals$seqname))
  if (isTRUE(check) & file.exists(outf)){
    ex_nms <- readLines(outf)[c(TRUE, FALSE, FALSE, FALSE)]
    if (nms %in% ex_nms){
      warning(sprintf("File %s already contains a read named %s", 
                      outf, nms))
    }
  }
  o <- file(outf, "a")
  seqname <- sprintf("@%s", nms)
  qualname <- sprintf("+%s", nms)
  writeLines(c(seqname, vals$seq, qualname, vals$quals) , o)
  close(o)
  return()
}
