#'Read a file in ab1 (Sanger) format and convert to fastq
#'
#'This is an R implementation of Wibowo Arindrarto's abifpy.py trimming module,
#'which itself implement's Richard Mott's trimming algorithm
#'See \url{https://github.com/bow/abifpy} for more details.
#'
#'Requires Bioconductor package SangerseqR
#'@param seqname name of sequence, to appear in fastq file
#'@param fname filename of sequence in ab1 format
#'@param outfname filename to append the fastq output to
#'@param trim should low quality bases be trimmed from the ends?  TRUE or FALSE
#'@param cutoff probability cutoff
#'@param min_seq_len minimum number of sequenced bases required in order to trim the read
#'@param offset phred offset for quality scores
#'@param recall Use sangerseqR to resolve the primary sequence if two sequences
#'are present.  May cause quality scores to be ignored. (Default: FALSE)
#'@author Helen Lindsay
#'@return None.  Sequences are appended to the outfname.
#'@examples
#'ab1_fname <- system.file("extdata", "IM2033.ab1", package = "CrispRVariants")
#'abifToFastq("IM2033", ab1_fname, "IM2033.fastq")
#'@export
abifToFastq <- function(seqname, fname, outfname, trim = TRUE, cutoff = 0.05,
                        min_seq_len = 20, offset = 33, recall = FALSE){
  sangerseqr <- requireNamespace("sangerseqR")
  stopifnot(isTRUE(sangerseqr))

  # Translation of python function
  abif <- sangerseqR::read.abif(fname)
  if (is.null(abif@data$PCON.2)){
    message(sprintf("failed on %s", seqname))
    return()
  }

  # Remove the extra character if it exists
  nucseq <- substring(abif@data$PBAS.2, 1, length(abif@data$PLOC.2))

  # sangerSeqR PCON.2 is a UTF8-encoded character vector in release, 
  # an integer vector in devel
  if (! typeof(abif@data$PCON.2) == "integer"){
    num_quals <- utf8ToInt(abif@data$PCON.2)[1:length(abif@data$PLOC.2)]
  } else {
    num_quals <- abif@data$PCON.2[1:length(abif@data$PLOC.2)]
  }
  
  if (isTRUE(recall)){
    recalled <- sangerseqR::makeBaseCalls(sangerseqR::sangerseq(abif))
    nucseq <- sangerseqR::primarySeq(recalled, string = TRUE)
    if (nchar(nucseq) != length(num_quals)){
      trim <- FALSE
      # Set all quality scores equal
      num_quals <- rep(60, nchar(nucseq)) 
      # 60 is compatible with all phred offsets, according to Wikipedia
      # Sanger scores do not usually exceed 60
      warning("Length of quality scores does not equal length of
              re-called base sequence, ignoring quality scores")
    }
  }
  
  if (trim == FALSE){
    writeFastq(outfname, list("seqname" = seqname, "seq" = nucseq,
                              "quals" = rawToChar(as.raw(num_quals + offset))))
    return()
  }

  trim_msg <- 'Sequence %s can not be trimmed because it is shorter than the trim
               segment size'
  if (nchar(nucseq) <= min_seq_len){
    warning(sprintf(trim_msg, seqname))
    return()
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

  # Additional check that there is enough sequence (not in abifpy):
  if (trim_finish - trim_start < min_seq_len -1){
    warning(sprintf(trim_msg, seqname))
    return()
  }

  writeFastq(outfname, list("seqname" = seqname,
              "seq" = substring(nucseq, trim_start, trim_finish),
              "quals" = rawToChar(as.raw(num_quals[trim_start:trim_finish]+offset))))
  return()
}
