#'@title Trims reads to a target region.
#'@description Trims aligned reads to a target region, optionally reverse
#'complementing the alignments.         
#'@param reads A GAlignments object, or a character
#'@param target A GRanges object specifying the range to narrow alignments to 
#'@param reverse.complement Should the alignments be oriented to match 
#'the strand of the target? (Default: TRUE)
#'@param collapse.pairs  If reads are paired, should pairs be collapsed?  (Default: FALSE) 
#'@param verbose Print progress and statistics (Default: FALSE)
#'@author Helen Lindsay
#'@rdname readsToTarget
#'
setGeneric("readsToTarget", function(reads, target, ...) {
  standardGeneric("readsToTarget")})


#'@param reads A GAlignments object
#'@rdname readsToTarget
setMethod("readsToTarget", signature("GAlignments", "GRanges"),
          function(reads, target, ..., reverse.complement = TRUE, collapse.pairs = FALSE, 
                   verbose = FALSE){
            
            
            
            if (reverse.complement & strand(target) == "*"){
             message(paste0("Target does not have a strand, but reverse.complement is TRUE.",
                            "Orienting reads to reference strand."))
            }else{
              rc <- rc.alns(strand(target), reverse.complement)
            }

            # MAKE SURE ALIGNMENTS ARE READ IN INCLUDING NAME AND FLAG IF REQUIRED
            
              #crun <- CrisprRun(bam, name)
            
            # After narrowing the alignments, collapse pairs
     
            
            
          })

#'
#'@param reads A character vector containing bam filenames to load
#'@rdname readsToTarget
setMethod("readsToTarget", signature("character", "GRanges"),
          function(reads, target, ..., reverse.complement = TRUE, collapse.pairs = FALSE, 
                   verbose = FALSE){  
            # lapply to character vector
          })



#'@title Internal crispRvariants function for deciding whether to reverse complement an alignment
#'@param target.strand Strand of the target region
#'@param reverse.complement Should the alignment be oriented to match the strand
#'@return A logical value indicating whether the narrowed alignment should be reverse complemented.
#'@author Helen Lindsay
rc.alns <- function(target.strand, reverse.complement){
  if (target.strand == "-" & reverse.complement = TRUE) return(TRUE)
  return(FALSE)
}


#'@title Internal crispRvariants function for collapsing pairs with concordant indels
#'@description Given a set of alignments to a target region, finds read pairs.
#'Compares insertion/deletion locations within pairs using the cigar string.  
#'Pairs with non-identical indels are excluded.  Pairs with identical indels are
#'collapsed to a single read, taking the consensus sequence of the pairs.
#'@param alns A GAlignments object.  We do not use GAlignmentPairs because amplicon-seq
#'can result in pairs in non-standard pairing orientation.  
#'Must include BAM flag, must not include unmapped reads.
#'@param use.consensus Should the consensus sequence be used if pairs have a mismatch?
#' (Default: TRUE)
#'@param keep.unpaired Should unpaired and chimeric reads be included?  (Default: TRUE)
#'@param verbose Report statistics on reads kept and excluded
#'@value The alignments, with non-concordant pairs removed and concordant pairs 
#'represented by a single read.
#'@author Helen Lindsay
collapse.pairs <- function(alns, use.consensus = TRUE, keep.unpaired = TRUE,
                           verbose = TRUE){  
  # 1 = 2^0 = paired flag
  # 2048 = 2^11 = supplementary alignment flag
  # As filtering for primary alignments occurs before finding duplicates, 
  # singletons are excluded
  is_primary <- !(bitwAnd(mcols(alns)$flag, 2048) & bitwAnd(mcols(alns)$flag, 1)) 
  pairs <- findChimeras(alns[is_primary])
  nms <- rle(names(alns)[is_primary][pairs]) 
  nms_codes <- rep(1:length(nms$lengths), nms$lengths)
  
  # If reads have the same insertions and deletions, they have identical cigar strings
  cig_runs <- rle(paste(cigar(alns)[is_primary][pairs], nms_codes, sep = "."))$lengths
  concordant <- rep(cig_runs, cig_runs) == rep(nms$lengths,nms$lengths)
  
  # Keep first alignment from all concordant pairs
  # Flag 64 = 2^6 = first alignment in pair
  is_pair <- which(is_primary)[pairs]
  is_first <- as.logical(bitwAnd(mcols(alns)$flag[is_pair], 64))
  keep <- is_pair[concordant & is_first]
  
  if (verbose){
    nunpaired <- length(alns) - length(is_pair)
    cc_t <- table(concordant)/2
    stats <- paste0("%s original alignments\n", 
              "  %s are not part of a primary alignment pair\n",
              "     (singletons and chimeras)\n",
              "  %s reads are paired \n",
              "    %s pairs have the same insertions/deletions\n",
              "    %s pairs have different insertions/deletions\n",
              "Keeping the first member of %s concordant read pairs\n")
    cat(sprintf(stats, length(alns), nunpaired, length(is_pair),
                cc_t[["TRUE"]], cc_t[["FALSE"]], cc_t[["TRUE"]]))
  }
  if (keep.unpaired){
    # Keep non-pairs, including non-primary and singletons
    keep <- c(keep, setdiff(c(1:length(alns)),is_pair))
    keep_alns <- alns[keep]
    if (verbose) cat(sprintf("Keeping %s unpaired reads\n", nunpaired))
  }
  
  if (use.consensus){  
    seq_runs <- rle(paste0(nms_codes, mcols(alns[pairs])$seq))$lengths
    same_seq <- rep(seq_runs, seq_runs) == rep(nms$lengths, nms$lengths)
 
    ncc_seqs <- mcols(alns[pairs][concordant & !same_seq])$seq
    consensus <- sapply(seq(1,length(ncc_seqs), by = 2), function(i){
      Biostrings::consensusString(ncc_seqs[i:(i+1)])
    })   
    # Overwrite the sequence of the non-concordant pairs.
    # The concordant alignments are at the start of keep
    #(concordant & ! same_seq)[seq(1,length(pairs), by = 2)]
    ncc_idxs <- cumsum(concordant & is_first)[concordant & is_first & !same_seq]
    mcols(keep_alns[ncc_idxs])$seq <- DNAStringSet(consensus)
    
    if (verbose){
      cat(sprintf("Using consensus for %s pairs with mismatches\n", 
                  length(consensus)))
    }
  }
  return(keep_alns)
}









  
  
  
  
initialize = function(bam, rc = FALSE, 
                      name = NULL, count_unmapped = TRUE, merge_chimeras = FALSE, 
                      tag_chimeras = TRUE, exclude_chimeras = FALSE, maxn = 100,
                      exclude_ranges = GRanges(), exclude_names = NA){


initialize = function(bam_fnames, ref, rc = FALSE,
                      short_cigars = FALSE, 
                      names = NULL, exclude_ranges = GRanges(), 
                      exclude_names = NA, 
                      merge_chimeras = TRUE, 
                      renumbered = FALSE, 
                      target_loc = NA, 
                      match_label = "no variant", 
                      upstream = 8,
                      downstream = 5, ...){
