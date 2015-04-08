#'Creates a text alignment from a set of cigar strings
#'
#'Creates a one-to-one text alignment of a set of cigar strings with respect
#'to the reference sequence by collapsing insertions and introducing gaps
#'across deletions. 
#'
#'When genomic coordinates for the alignment start and the target region 
#'are provided, aligned sequences are cropped to the target region 

#'@author Helen Lindsay
#'@param cigar A list of cigar strings to align
#'@param dnaseq The set of sequences corresponding to the cigars, as DNAStrings
#'@param del_char The character to represent deleted bases. Default "-" 
#'@param aln_start Genomic start locations of aligned sequences. Should be
#'used in conjunction with target_start and target_end.
#'@param target_start Genomic start of the region to be returned. 
#'@param target_end Genomic end of the region to be returned.
#'@return The sequences with insertions collapsed and deletions padded
#'
seqsToAln <- function(cigar, dnaseq, del_char = "-", aln_start = NULL, target_start = NULL, 
                      target_end = NULL){
  # TO DO - CHECK THAT THIS IS VECTORISED
  
  # Remove insertion sequences
  wrt_ref <- cigarRangesAlongReferenceSpace(cigar)[[1]]
  wrt_read <- cigarRangesAlongQuerySpace(cigar)[[1]]
  ops <- explodeCigarOps(cigar)[[1]]
  segs <- as.character(Views(dnaseq, wrt_read))
  segs[which(ops == "I")] <- ""
  for (j in which(ops == "D")){
    segs[j] <- paste0(rep(del_char, width(wrt_ref[j])), collapse = "")
  }
  result <- paste0(segs, collapse = "")
  
  if (! is.null(aln_start) & ! is.null(target_start) & ! is.null(target_end) ){
    trim_start <- target_start - (aln_start - 1)
    
    if ( trim_start < 0){
      stop("dnaseq to be trimmed must start before the target location")
    }
    
    trim_end <- trim_start + (target_end - target_start)
    if (trim_end < length(result)){
      stop("dnaseq is not long enough to trim to the target region")
    }
    result <- subseq(result, trim_start,trim_end)
  }
  result
}