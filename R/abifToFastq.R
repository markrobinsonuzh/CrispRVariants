#'Read a file in ab1 (Sanger) format and convert to fastq
#'
#'This is an R implementation of Wibowo Arindrarto's abifpy.py trimming module,
#'which itself implement's Richard Mott's trimming algorithm
#'See \link{https://github.com/bow/abifpy} for more details.
#'
#'Requires Bioconductor package SangerseqR
#'@param seqname name of sequence, to appear in fastq file
#'@param fname filename of sequence in ab1 format
#'@param outfname filename to append the fastq output to
#'@param trim should low quality bases be trimmed from the ends?  TRUE or FALSE
#'@param cutoff probability cutoff
#'@param min_seq_len minimum number of sequenced bases required in order to trim the read
#'@param offset phred offset for quality scores
#'@author Helen Lindsay
#'@export
abifToFastq <- function(seqname, fname, outfname, trim = TRUE, cutoff = 0.05, 
                        min_seq_len = 20, offset = 33){
  # TO DO? check min_seq_len after trimming also?
  
  sangerseqr <- require(sangerseqR)
  stopifnot(sangerseqr == TRUE)
  
  # Translation of python function
  abif <- sangerseqR::read.abif(fname)
  if (is.null(abif@data$PCON.2)){
    print(sprintf("failed on %s", seqname))
    return()
  }
  
  # Remove the extra character if it exists
  nucseq <- substring(abif@data$PBAS.2, 1, length(abif@data$PLOC.2))
  num_quals <- utf8ToInt(abif@data$PCON.2)[1:length(abif@data$PLOC.2)] 
  
  if (trim == FALSE){
    writeFastq(outfname, list("seq" = nucseq, "quals" = rawToChar(as.raw(num_quals + 33))))
    return()
  }
  
  if (nchar(nucseq) <= min_seq_len){
    stop('Sequence can not be trimmed because it is shorter than the trim segment size')
  } 
  scores = cutoff - 10^(num_quals / -10)
  running_sum <- rep(0, length(scores) + 1)
  
  # Note running_sum counts from zero, scores from 1 
  for (i in 1:length(scores)){
    num <- scores[i] + running_sum[i]
    running_sum[i+1] <- ifelse(num < 0, 0, num)
  }
  
  trim_start <- min(which(running_sum > 0)) - 1
  trim_finish <- which.max(running_sum) - 2 
  # -1 for running_sum offset, -1 because python doesn't include ends
  
  writeFastq(outfname, list("seqname" = seqname, "seq" = substring(nucseq, trim_start, trim_finish),
                            "quals" = rawToChar(as.raw(num_quals[trim_start:trim_finish]+offset))))
  return()
}
