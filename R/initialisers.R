#'@title Trims reads to a target region.
#'@description Trims aligned reads to a target region, optionally reverse
#'complementing the alignments.         
#'@param reads A GAlignments object, or a character vector of the filenames
#'@param target A GRanges object specifying the range to narrow alignments to 
#'@param reverse.complement Should the alignments be oriented to match 
#'the strand of the target? (Default: TRUE)
#'@param collapse.pairs  If reads are paired, should pairs be collapsed?  (Default: FALSE) 
#'@param keep.unpaired For paired reads, should unpaired alignments be kept? (Default: TRUE)
#'@param verbose Print progress and statistics (Default: FALSE)
#'@author Helen Lindsay
#'@rdname readsToTarget
#'@export
setGeneric("readsToTarget", function(reads, target, ...) {
  standardGeneric("readsToTarget")})


#'@param name An experiment name for the reads.  (Default: NULL)  
#'@rdname readsToTarget
setMethod("readsToTarget", signature("GAlignments", "GRanges"),
          function(reads, target, ..., reverse.complement = TRUE, 
                   collapse.pairs = FALSE, verbose = FALSE, name = NULL){
          
            if (length(target) > 1){
              stop("readsToTarget accepts a single target range")
            }
            if (collapse.pairs = TRUE){
              if (is.null(names(reads)) |  !"flag" in names(mcols(reads))){
                stop("Reads must include names and bam flags for collapsing pairs")
              }
            }  
            # Filter out reads that don't span the target region 
            # Not using findOverlaps because reads may be paired, i.e. names nonunique
            bam <- reads[start(reads) <= start(target) & end(reads) >= end(target) & 
                         seqnames(reads) == as.character(seqnames(target))]
            
            if (verbose){
              cat(sprintf("%s of %s reads are on target\n", length(bam), length(reads)))
            }
            
            if (reverse.complement & strand(target) == "*"){
             message(paste0("Target does not have a strand, but reverse.complement is TRUE.",
                            "Orienting reads to reference strand."))
            }else{
              rc <- rc.alns(strand(target), reverse.complement)
            }
            
            # narrow aligned reads
            result <- narrowAlignments(bam, target, reverse.complement = rc, verbose)
            bam <- result$alignments
            
            # Collapse pairs of narrowed reads
            result <- collapsePairs(bam, genome.ranges = result$genome.ranges)
            bam <- result$alignments
            genome.ranges <- result$genome.ranges
              #crun <- CrisprRun(bam, name)      
            
          })

#'@param names Experiment names for each bam file.  If not supplied, filenames are used.
#'@param chimeras Flag to determine how chimeric reads are treated.  One of 
#'"ignore", "exclude", and "merge".  Default "ignore", "merge" not implemented yet
#'@rdname readsToTarget
setMethod("readsToTarget", signature("character", "GRanges"),
          function(reads, target, ..., reverse.complement = TRUE, 
                   exclude.ranges = GRanges(), exclude.names = NA,
                   chimeras = c("ignore","exclude","merge")
                   collapsePairs = FALSE, names = NULL, verbose = FALSE){  
            
            alns <- lapply(reads, readTargetBam, exclude.ranges = exclude.ranges,
                           exclude.names = exclude.names, chimeras = chimeras,
                           verbose = verbose)
            #crispr.runs <- lapply(alns, readsToTarget) 
          })


#'@title Internal crispRvariants function for reading and filtering a bam file
#'@description Includes options for excluding reads either by name or range.
#'The latter is useful if chimeras are excluded.  Reads are excluded before
#'chimeras are detected, thus a chimeric read consisting of two sections, one of 
#'which overlaps an excluded region, will not be considered chimeric.
#'Chimeric reads can be ignored, excluded, which means that all sections of a 
#'chimeric read will be removed, or merged, which means that chimeras will be
#'collapsed into a single read where possible. (Not implemented yet)
#'If chimeras = "merge", chimeric reads are merged if all segments
# are from the same chromosome, do not overlap, and are aligned to the same strand.
# It is assumed that sequences with two alignments are chimeras, not alternate mappings
#'@param file The name of a bam file to read in
#'@param target A GRanges object containing a single target range
#'@param exclude.ranges A GRanges object of regions that should not be counted,
#'e.g. primer or cloning vector sequences that have a match in the genome 
#'@param exclude.names A vector of read names to exclude.  
#'@param chimeras Flag to determine how chimeric reads are treated.  One of 
#'"ignore", "exclude", and "merge".  Default "ignore".  
readTargetBam <- function(file, target, exclude.ranges = GRanges(), 
                          exclude.names = NA, 
                          chimeras = c("ignore","exclude","merge")){
  
  ch.action <- match.arg(chimeras, c("ignore","exclude","merge"))
  if (ch.action == "ignore"){
    # If chimeras are not to be excluded or merged, 
    # we only need to read in reads overlapping the target region
    param <- ScanBamParam(what = c("seq", "flag"), which = target)
  } else {
    # In this case, must read in the entire bam to be sure of finding chimeric reads
    param <- ScanBamParam(what = c("seq", "flag"))
  }
  bam <- GenomicAlignments::readGAlignments(file, param = param, use.names = TRUE)
  #Exclude reads by name or range
  temp <- excludeFromBam(bam, exclude.ranges, exclude.names)    
  
  if (verbose == TRUE){
    original <- length(bam)
    cat(sprintf("Read %s alignments, excluded %s\n", original, original - length(temp)))
  }
  bam <- temp 
  
  if (length(bam) == 0 | ch.action == "ignore") return(bam)
  
  chimera_idxs <- findChimeras(bam)
  
  if (chimeras == "exclude"){
    if( length(chimera_idxs) >= 2){
      bam <- bam[-chimera_idxs]
    }
    if (verbose == TRUE){
      cat(sprintf("%s reads after filtering chimeras\n", length(bam)))
    }
    return(bam)
  }
  if (chimeras == "merge"){
    cat("Merging chimeras not implemented yet, ignoring chimeras") 
    return(bam)
  }
}


#'@title Internal crispRvariants function for deciding whether to reverse 
#'complement aligned reads
#'@param target.strand Strand of the target region
#'@param reverse.complement Should the alignment be oriented to match the strand
#'@return A logical value indicating whether the narrowed alignment should be reverse complemented.
#'@author Helen Lindsay
rc.alns <- function(target.strand, reverse.complement){
  if (target.strand == "-" & reverse.complement == TRUE) return(TRUE)
  return(FALSE)
}


#'@title Narrow a set of aligned reads to a target region
#'@description Aligned reads are narrowed to the target region.  In
#'the case of reads with deletions spanning the boundaries of the target,
#'reads are narrowed to the next aligned base outside of the target
#'
#'Note that alignments and cigars are reversed if reverse.complement = TRUE
# but the start is still genomic. i.e. w.r.t. the reference strand

#'@param alns A GAlignments object including a metadata column "seq" 
#'containing the sequence
#'@param target A GRanges object
#'@param reverse.complement Should the aligned reads be reverse complemented?
#'@param verbose (Default: FALSE)
#'@author Helen Lindsay
#'@rdname narrowAlignments
setGeneric("narrowAlignments", function(alns, target, ...) {
  standardGeneric("narrowAlignments")})

#'@return A list of "alignments", the narrowed alignments (GAlignments), 
#'and "genome.ranges", the genomic locations of the cigar operations with
#'respect to the reference strand.
#'If reverse.complement is true, "alignments" is not a proper GAlignments object.
#'The start is displayed with respect to the reference strand, but the cigar
#'and sequence are displayed with respect to the negative strand.  
#'Use with caution outside of crispRvariants!
#'@rdname narrowAlignments 
setMethod("narrowAlignments", signature("GAlignments", "GRanges"),
          function(alns, target, ..., reverse.complement, verbose = FALSE){  
            # Narrowing example:
            # 3-4-5-6-7-8-9-10 Read
            #     5-6-7-8      Target sequence
            # target_start 5 - (read_start 3 - 1) = index 3
            # target_end 8 - target_start 5 + cigstart 3 = index 6
            ref.ranges <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar(alns))
            query.ranges <- GenomicAlignments::cigarRangesAlongQuerySpace(cigar(alns))
            cig.ops <- GenomicAlignments::explodeCigarOps(cigar(alns))
            cig.ops <- IRanges::CharacterList(cig.ops)
            if (verbose == TRUE) cat("finding deletions \n")
            locs <- findDeletions(start(target), end(target), alns, ref.ranges, cig.ops)
            if (verbose == TRUE) cat("narrowing alignments\n")
            
            # Record how many bases before first aligned read    
            clip_starts <- rep(0, length(alns))
            is_clipped <- grepl("^[0-9]+[HS]", cigar(alns))
            clip_starts[is_clipped] <- unlist(lapply(width(query.ranges[is_clipped]), '[[', 1))
            
            temp <- GenomicAlignments::cigarNarrow(cigar(alns), locs$starts, locs$ends)  
            new_starts <-  attr(temp, "rshift") + 1 + clip_starts 
            # + 1 as rshift is number removed not starting point
            # new_starts are read offsets wrt current genomic starts.
            
            # Account for insertions and deletions prior to the target site
            # Need to adjust for deletions as new_starts derived wrt reference
            prior <- start(ref.ranges) <= locs$start
            prior_ins <- prior & cig.ops == "I"
            prior_del <- prior & cig.ops == "D" 
            new_starts <- new_starts + sum(width(query.ranges)[prior_ins])
            new_starts <- new_starts - sum(width(ref.ranges)[prior_del])
          
            cigs <- as.character(temp)
            query.ranges <- GenomicAlignments::cigarRangesAlongQuerySpace(cigs)
            shift_starts <- start(alns) + attr(temp, "rshift") -1
            gr <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigs)
            genome_ranges <- shift(gr, shift_starts)
            
            seq_lens <- sum(width(query.ranges))  
            seqs <- subseq(mcols(alns)$seq, start = new_starts, width = seq_lens)
            genome_start <- start(alns) + attr(temp, "rshift")
            
            if (reverse.complement == TRUE){
              if (verbose == TRUE) cat("reversing narrowed alignments\n")
              cigs <- unname(sapply(cigs, reverseCigar))
              query.ranges <- GenomicAlignments::cigarRangesAlongQuerySpace(cigs)
              
              # Note: even if strand is -ve, it is displayed wrt reference strand in bam
              seqs <- Biostrings::reverseComplement(seqs)
              genome_ranges <- IRanges::IRangesList(lapply(.self$genome_ranges, rev))
              
              # Add difference between target start and actual start, subtract difference at end
              lshift <- target_start - (start(alns) - 1) - locs$starts
              rshift <- target_end - target_start - (locs$ends - locs$starts - lshift)           
              genome_start <- as.integer(genome_start + lshift + rshift)
            }
            
            ga_params <- list(seqnames = seqnames(alns), pos = genome_start, cigar = cigs, 
                              names = names(alns), strand = strand(alns), 
                              seqlengths = seqlengths(alns), seq = seqs)
            
            mcols <- as.list(mcols(alns))
            mcols$seq <- NULL
            ga_params <- c(ga_params, mcols)
            alns <- do.call(GAlignments, ga_params) 
            
            return(list(alignments = alns, "genome.ranges" = genome_ranges))
          })


#'@title Find location of deletions that overlap a target region.  
#'@description   For reads with a deletion spanning one or both
#' ends of the target location, narrow the alignment to encompass the deletion 
#'Deletions may be coded as either "D" or "N" (splice junction), 
# depending upon the mapping software. 
#'@param target.start The target start location
#'@param target.end The target end location
#'@param alns GAlignments, only containing alignments spanning the target region
#'@param ref.ranges IRangesList created with GenomicAlignments::cigarRangesAlongReferenceSpace
#'@param cig.ops Cigar operations, created with GenomicAlignments::explodeCigarOps
#'@param del.chars Characters that may represent deletions.  Default: c("D", "N") 
#'@return A list of the start and end locations w.r.t. the reads for passing to 
#'GenomicAlignments::CigarNarrow
findDeletions <- function(target.start, target.end, alns, ref.ranges, 
                          cig.ops, del.chars = c("D", "N")){
  
  idxs <- rep(1:length(cig.ops), lapply(cig.ops, length))
  genomic <- shift(ref.ranges, start(alns)-1) 
  on_target <- unlist(start(genomic) <= target.end & end(genomic) >= target.start)
  codes <- paste0(idxs, on_target)
  s_del <- !duplicated(codes) & on_target & unlist(cig.ops) %in% del.chars
  e_del <- rev(!duplicated(rev(codes))) & on_target & unlist(cig.ops) %in% del.chars
  
  # Get the ranges for narrowing in read coordinates
  result_s <- target.start - (start(alns) - 1)
  result_e <- target.end - target.start  + result_s  
  result_s[idxs[s_del]] <- unlist(start(ref.ranges))[s_del] - 1
  result_e[idxs[e_del]] <- unlist(end(ref.ranges))[e_del] + 1
  
  return(list(starts = result_s, ends = result_e)) 
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
#'Setting this to be TRUE makes this function much slower (Default: TRUE)
#'@param keep.unpaired Should unpaired and chimeric reads be included?  (Default: TRUE)
#'@param verbose Report statistics on reads kept and excluded
#'@param ... Additional items with the same length as alns, 
#'that should be filtered to match alns.
#'@return The alignments, with non-concordant pairs removed and concordant pairs 
#'represented by a single read.
#'@author Helen Lindsay
collapsePairs <- function(alns, use.consensus = TRUE, keep.unpaired = TRUE,
                          verbose = TRUE, ...){  
  dots <- list(...)
  print(dots)
  if (! unique(sapply(dots, length)) == length(alns)){
    stop("Each ... argument supplied must have the same length as the alignments")
  }
  
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
    filtered.dots <- lapply(dots, function(x) x[keep])    
  }
  return(c(list("alignments" = keep_alns), filtered.dots))
}

