#' @title Creates a text alignment from a set of cigar strings
#'@description Creates a one-to-one text alignment of a set of cigar strings with respect
#'to the reference sequence by collapsing insertions and introducing gaps
#'across deletions.
#'
#'When genomic coordinates for the alignment start and the target region
#'are provided, aligned sequences are cropped to the target region
#'@author Helen Lindsay
#'@param cigar A list of cigar strings to align
#'@param dnaseq The set of sequences corresponding to the cigars, as Biostrings::DNAStrings
#'@param target The target region to return, as GRanges.  Sequences overlapping
#'the target region are trimmed to exactly match it.
#'@param del_char The character to represent deleted bases. Default "-"
#'@param aln_start Genomic start locations of aligned sequences. Should be
#'used in conjunction with target_start and target_end.
#'@return The sequences with insertions collapsed and deletions padded
#'@rdname seqsToAln
seqsToAln <- function(cigar, dnaseq, target, del_char = "-", aln_start = NULL){
  # Additional trimming is necessary because deletion operations may overhang
  # boundaries of target region
  target_start <- start(target)
  target_end <- end(target)
  strand <- as.character(strand(target))
  if (strand == "*") strand <- "+"
  
  result <- GenomicAlignments::sequenceLayer(dnaseq, cigar, D.letter = del_char,
                                          N.letter = "N")
  sq_len <- GenomicAlignments::width(result)

  if (! is.null(aln_start) & ! is.null(target_start) & ! is.null(target_end) ){
    trim_start <- target_start - (aln_start - 1)

    if (any(trim_start < 0)){
      stop("dnaseq to be trimmed must start before the target location")
    }

    trim_end <- trim_start + (target_end - target_start)
    if (any(trim_end > sq_len)){
      stop("dnaseq is not long enough to trim to the target region")
    }
    result <- subseq(result, trim_start,trim_end)
  }
  if (strand == "-"){
    result <- reverseComplement(result)
  }
  result <- as.character(result)
  result
}