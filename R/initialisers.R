#'@title Trims reads to a target region.
#'@description Trims aligned reads to one or several target regions, 
#'optionally reverse complementing the alignments.         
#'@param reads A GAlignments object, or a character vector of the filenames
#'@param target A GRanges object specifying the range to narrow alignments to 
#'@author Helen Lindsay
#'@rdname readsToTarget
#'
#'@import BiocParallel
#'@import Biostrings
#'@import ggplot2
#'@import grid
#'@import gridExtra
#'@import GenomicAlignments
#'@import GenomicRanges
#'@import IRanges
#'@import methods
#'@import Rsamtools
#'@importFrom reshape2 melt
#'@importFrom AnnotationDbi select
#'@export
setGeneric("readsToTarget", function(reads, target, ...) {
  standardGeneric("readsToTarget")})


#'@param name An experiment name for the reads.  (Default: NULL)  
#'@param reverse.complement Should the alignments be oriented to match 
#'the strand of the target? (Default: TRUE)
#'@param collapse.pairs  If reads are paired, should pairs be collapsed?  (Default: FALSE) 
#'Note: only collapses primary alignments, and assumes that there is only one primary
#'alignment per read.  May fail with blat alignments converted to bam. 
#'@param use.consensus Take the consensus sequence for non-matching pairs? If FALSE,
#'the sequence of the first read is used.  Can be very slow. (Default: FALSE)
#'@param store.chimeras Should chimeric reads be stored?  (Default: FALSE)
#'@param verbose Print progress and statistics (Default: FALSE)
#'@return (signature("GAlignments", "GRanges")) A \code{\link{CrisprRun}} object 
#'@examples
#'# Load the metadata table
#'md_fname <- system.file("extdata", "gol_F1_metadata_small.txt", package = "crispRvariants")
#'md <- read.table(md_fname, sep = "\t", stringsAsFactors = FALSE)
#'
#'# Get bam filenames and their full paths
#'bam_fnames <- sapply(md$bam.filename, function(fn){
#'  system.file("extdata", fn, package = "crispRvariants")})
#'
#'reference <- DNAString("GGTCTCTCGCAGGATGTTGCTGG")
#'gd <- GRanges("18", IRanges(4647377, 4647399), strand = "+")
#'
#'crispr_set <- readsToTarget(bam_fnames, target = gd, reference = reference,
#'                            names = md$experiment.name, target.loc = 17)
#'
#'@rdname readsToTarget
setMethod("readsToTarget", signature("GAlignments", "GRanges"),
          function(reads, target, ..., reverse.complement = TRUE, 
                   collapse.pairs = FALSE, use.consensus = FALSE, 
                   store.chimeras = FALSE, verbose = FALSE, name = NULL){
            
            if (length(target) > 1){
              stop("readsToTarget accepts a single target range")
            }
            if (collapse.pairs == TRUE){
              if (is.null(names(reads)) |  ! ("flag" %in% names(mcols(reads))) ){
                stop("Reads must include names and bam flags for collapsing pairs")
              }
            }  
            # Filter out reads that don't span the target region 
            # Not using findOverlaps because reads may be paired, i.e. names nonunique
            
            bam <- reads[start(reads) <= start(target) & end(reads) >= end(target) & 
                         seqnames(reads) == as.character(seqnames(target))]
            
            if (store.chimeras == TRUE){
              # ASSUME THAT TARGET CHIMERAS HAVE BEEN MERGED AT THIS POINT
              # find chimeric alignments for the reads that span the target 
              # this step ensures that off-target reads are not kept
              ch <- reads[findChimeras(reads[names(reads) %in% names(bam)])]
              # which chimeras are not already in the bam
              is_ch <- ! granges(ch) %in% granges(bam)
              chimeras <- ch[is_ch]
            } else {
              chimeras <- GAlignments()
            }
              
            if (verbose){
              cat(sprintf("%s of %s reads span the target range\n", length(bam), length(reads)))
            }
            if (length(bam) == 0) return(NULL)
            
            if (reverse.complement & as.character(strand(target)) == "*"){
             message(paste0("Target does not have a strand, but reverse.complement is TRUE.",
                            "Orienting reads to reference strand."))
              rc = FALSE
            }else{
              rc <- rcAlns(as.character(strand(target)), reverse.complement)
            }
            
            # narrow aligned reads
            result <- narrowAlignments(bam, target, reverse.complement = rc, 
                                       verbose = verbose)
            bam <- result$alignments
            
            # Collapse pairs of narrowed reads
            if (collapse.pairs == TRUE){
              result <- collapsePairs(bam, genome.ranges = result$genome.ranges, 
                                      use.consensus = use.consensus, verbose = verbose)
            }
            if (is.null(result)) return(NULL)
            
            bam <- result$alignments
            genome.ranges <- result$genome.ranges
            
            crun <- CrisprRun(bam, target, genome.ranges, rc = rc, name = name, 
                              chimeras = chimeras, verbose = verbose)      
            return(crun)
          })


#'@param names Experiment names for each bam file.  If not supplied, filenames are used.
#'@param chimeras Flag to determine how chimeric reads are treated.  One of 
#'"ignore", "exclude", and "merge".  Default "ignore", "merge" not implemented yet
#'@param reference The reference sequence
#'@param exclude.ranges Ranges to exclude from consideration, e.g. homologous to a pcr primer.
#'@param exclude.names Alignment names to exclude
#'@param ... Extra arguments for initialising CrisprSet
#'@return (signature("character", "GRanges")) A \code{\link{CrisprSet}} object
#'@rdname readsToTarget
setMethod("readsToTarget", signature("character", "GRanges"),
          function(reads, target, ..., reference, reverse.complement = TRUE, 
                   exclude.ranges = GRanges(), exclude.names = NA,
                   chimeras = c("ignore","exclude","merge"),
                   collapse.pairs = FALSE, use.consensus = TRUE, 
                   names = NULL, verbose = FALSE){  
            
            # Make sure the reference sequence consists of one sequence
            # Coerce to "DNA string if necessary"
            if (class(reference) == "character" | class(reference) == "DNAStringSet"){
              if (length(reference) > 1){
                stop("Reference should be a single sequence, as a DNAString")
              }
              reference <- as(reference[1], "DNAString")
            }
            
            alns <- lapply(reads, readTargetBam, target = target, 
                           exclude.ranges = exclude.ranges,
                           exclude.names = exclude.names, chimeras = chimeras,
                           verbose = verbose)
            
            # If names are not specified, set them to the filenames
            if (is.null(names)){
              names <- reads
            }
            names <- as.character(names)
            cset <- alnsToCrisprSet(alns, reference, target, reverse.complement, 
                                    collapse.pairs, names, use.consensus, 
                                    verbose, store.chimeras, ...)
            return(cset)
          })


#'@export
#'@rdname readsToTarget
setGeneric("readsToTargets", function(reads, targets, ...) {
  standardGeneric("readsToTargets")})

#'@param targets A set of targets to narrow reads to
#'@param references A set of reference sequences matching the targets
#'@param primer.ranges A set of GRanges, corresponding to the targets.
#'Read lengths are typically greater than target regions, and it can 
#'be that reads span multiple targets.  If primer.ranges are available,
#'they can be used to assign such reads to the correct target.  
#'@param ignore.strand Should strand be considered when finding overlaps?
#'(See \code{\link[GenomicAlignments]{findOverlaps}} )
#'@param bpparam A BiocParallel parameter for parallelising across reads.  
#'Default: no parallelisation.  (See \code{\link[BiocParallel]{bpparam}})
#'@rdname readsToTarget
setMethod("readsToTargets", signature("character", "GRanges"),
          function(reads, targets, ..., references, primer.ranges = NULL, 
                   reverse.complement = TRUE, collapse.pairs = FALSE, 
                   use.consensus = TRUE, ignore.strand = TRUE, 
                   names = NULL, bpparam = BiocParallel::SerialParam(), 
                   verbose = FALSE){
          
            # ACCOUNT FOR CHIMERIC READS OR NOT?
            
            if (! is.null(primer.ranges)){
              if (! (length(primer.ranges) == length(targets))){
                stop("primer.ranges should contain one range per target")
              }
            }
            if (is.null(names)){
              names <- reads
            }
            
            param <- Rsamtools::ScanBamParam(what = c("seq", "flag"))
            args <- list(...)
                      
            bamsByPCR <- bplapply(seq_along(reads), function(i){
              if (verbose) cat(sprintf("Loading alignments for %s\n\n", names[i]))
              
              bam <- GenomicAlignments::readGAlignments(reads[i], 
                                              param = param, use.names = TRUE)
              if (length(bam) == 0){
                if (verbose) cat("No reads in alignment\n")
                return(NULL)
              } 
         
              # If primer.ranges are provided, match reads to primers
              # If not, match reads to targets 
              if (! is.null(primer.ranges)){
                hits <- readsByPCRPrimer(bam, primer.ranges, verbose = verbose)
                #if (is.null(hits)) 
                bamByPCR <- split(bam[queryHits(hits)], subjectHits(hits))
              } else{
                hits <- findOverlaps(targets, bam, type = "within")
                duplicates <- (duplicated(subjectHits(hits)) | 
                                duplicated(subjectHits(hits), fromLast = TRUE))                
                if (verbose){
                  msg <- paste0("%s (%.2f%%) reads of %s overlap a target\n",
                                "  %s (%.2f%%) of these overlapping multiple targets removed\n",
                                "  %s (%.2f%%) reads mapped to a single target\n\n")
                  rhits <- length(unique(subjectHits(hits)))
                  bl <- length(bam)
                  ndups <- sum(duplicates)
                  nndups <- sum(!duplicates)
                  cat(sprintf(msg, rhits, rhits/bl*100, bl, ndups, ndups/rhits*100, 
                              nndups, nndups/bl*100))
                } 
                hits <- hits[!duplicates]
                bamByPCR <- split(bam[subjectHits(hits)], queryHits(hits))
              }
              bamByPCR
            }, BPPARAM = bpparam)
            
            #bamsByPCR <- bamsByPCR[ !(is.null)
            
            # bamsByPCR is separated by sample
            # rearrange to separate by target, initialise CrisprSets
            # optionally perform functions with these
            bbpcr_nms <- lapply(bamsByPCR, names)
            bbpcr_idx <- rep(1:length(bamsByPCR), lapply(bbpcr_nms, length))
            bbpcr <- do.call(c, unlist(bamsByPCR, use.names = FALSE))
            tgts <- unlist(bbpcr_nms, use.names = FALSE)
            unq_tgts <- unique(tgts)
            result <- bplapply(unq_tgts, function(tgt){
              if (verbose == TRUE) cat(sprintf("\n\nWorking on target %s\n", tgt))
              # which bams include this target
              idxs <- which(tgts == tgt)
              bams <- bbpcr[idxs]
              names(bams) <- names[bbpcr_idx[idxs]]
              target <- targets[as.numeric(tgt)]
              reference <- references[[as.numeric(tgt)]]
              cset <- alnsToCrisprSet(as.list(bams), reference, target, reverse.complement, 
                                      collapse.pairs, names(bams), use.consensus, 
                                      verbose, ...)
            }, BPPARAM = bpparam)
            
            if (! is.null(names(targets))){
              names(result) <- names(targets)[as.numeric(unq_tgts)]
            }
            result <- result[! sapply(result, is.null)]
            return(result)
          })



alnsToCrisprSet <- function(alns, reference, target, reverse.complement,
                            collapse.pairs, names, use.consensus, 
                            verbose, store.chimeras = FALSE, ...){
  print(sprintf("Processing %s samples", length(alns)))
  print(sprintf("verbose = %s", verbose))
  crispr.runs <- lapply(seq_along(alns), function(i){
    crun <- readsToTarget(alns[[i]], target = target, 
                reverse.complement = reverse.complement,
                collapse.pairs = collapse.pairs, use.consensus = use.consensus,
                verbose = verbose, name = names[i])
    crun
  })
  
  to_rm <- sapply(crispr.runs, is.null)
  if (any(to_rm)){
    if (verbose){
      rm_nms <- paste0(names[to_rm], collapse = ",", sep = "\n")
      cat(sprintf("Excluding samples that have no on target reads:\n%s",
                  rm_nms))
    }
    crispr.runs <- crispr.runs[!to_rm]
    names <- names[!to_rm]
    if (length(crispr.runs) == 0) {
      warning("Could not narrow reads to target, no samples have on-target alignments")
      return()
    }
  }
  
  rc <- rcAlns(as.character(strand(target)),reverse.complement)
  cset <- CrisprSet(crispr.runs, reference, target, rc = rc,
                    verbose = verbose, names = names, ...)
  return(cset)
}


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
#'@param verbose Print stats about number of alignments read and filtered.  (Default: FALSE)
#'@return A GenomicAlignments::GAlignment obj
readTargetBam <- function(file, target, exclude.ranges = GRanges(), 
                          exclude.names = NA,
                          chimeras = c("ignore","exclude","merge"),
                          verbose = FALSE){
  
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
rcAlns <- function(target.strand, reverse.complement){
  if (target.strand == "-" & reverse.complement == TRUE) return(TRUE)
  return(FALSE)
}


#'@title Narrow a set of aligned reads to a target region
#'@description Aligned reads are narrowed to the target region.  In
#'the case of reads with deletions spanning the boundaries of the target,
#'reads are narrowed to the next aligned base outside of the target
#'Note that alignments and cigars are reversed if reverse.complement = TRUE
# but the start is still genomic. i.e. w.r.t. the reference strand
#'@param alns A GAlignments object including a metadata column "seq" 
#'containing the sequence
#'@param target A GRanges object
#'@param reverse.complement Should the aligned reads be reverse complemented?
#'@param verbose (Default: FALSE)
#'@param ... additional arguments
#'@author Helen Lindsay
#'@rdname narrowAlignments
#'@export
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
            # target_end 8 - target_start 5 + cigstart 3 = index 
            
            target_start <- start(target)
            target_end <- end(target)
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
            
            reverse_ranges <- function(gr){
              # Quicker version of IRanges::IRangesList(lapply(list, rev))
              all_r <- rev(unlist(gr))
              sp <- rev(rep(1:length(gr), lapply(gr, length)))
              # Note: splitting works in numeric/factor order
              result <- as(split(all_r, sp), "IRangesList")
              names(result) <- NULL
              return(result)
            }
            
            if (reverse.complement == TRUE){
              if (verbose == TRUE) cat("reversing narrowed alignments\n")
              cigs <- unname(sapply(cigs, reverseCigar))
              query.ranges <- GenomicAlignments::cigarRangesAlongQuerySpace(cigs)
              
              # Note: even if strand is -ve, it is displayed wrt reference strand in bam
              seqs <- Biostrings::reverseComplement(seqs)
              genome_ranges <- reverse_ranges(genome_ranges)
                            
              # Add difference between target start and actual start, subtract difference at end
              lshift <- target_start - (start(alns) - 1) - locs$starts
              rshift <- target_end - target_start - (locs$ends - locs$starts - lshift)           
              genome_start <- as.integer(genome_start + lshift + rshift)
              if (verbose == TRUE) cat("alignments reversed\n")
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
  # find reads with a deletion spanning the target start or end
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
  
  if (! unique(sapply(dots, length)) == length(alns)){
    stop("Each ... argument supplied must have the same length as the alignments")
  }
  
  # 1 = 2^0 = paired flag
  # 2048 = 2^11 = supplementary alignment flag
  is_primary <- !(bitwAnd(mcols(alns)$flag, 2048) & bitwAnd(mcols(alns)$flag, 1)) 
  pairs <- findChimeras(alns[is_primary]) # This just matches read names
  
  # If there are no pairs, no need to do anything further 
  if (length(pairs) == 0){
    if (keep.unpaired == TRUE){
      return(c(list("alignments" = alns), dots))
    } else {
      return(NULL)
    }
  }
  # Pairs are primary alignments with the same name
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
    cc_true <- sum(concordant)/2
    cc_false <- sum(!concordant)/2
    stats <- paste0("\nCollapsing paired alignments:\n",
              "%s original alignments\n", 
              "  %s are not part of a primary alignment pair\n",
              "     (singletons and chimeras)\n",
              "  %s reads are paired \n",
              "    %s pairs have the same insertions/deletions\n",
              "    %s pairs have different insertions/deletions\n",
              "Keeping the first member of %s concordant read pairs\n")
    cat(sprintf(stats, length(alns), nunpaired, length(is_pair),
                cc_true, cc_false, cc_true))
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
    if (verbose){
      cat(sprintf("Finding consensus for %s pairs with mismatches\n", 
                  length(ncc_seqs)/2))
    }
    if (length(ncc_seqs) >= 2){
      consensus <- sapply(seq(1,length(ncc_seqs), by = 2), function(i){
      Biostrings::consensusString(ncc_seqs[i:(i+1)])
      })
      # Overwrite the sequence of the non-concordant pairs.
      # The concordant alignments are at the start of keep
      ncc_idxs <- cumsum(concordant & is_first)[concordant & is_first & !same_seq]
      mcols(keep_alns[ncc_idxs])$seq <- Biostrings::DNAStringSet(consensus)
    }
  }  
  if (length(keep) == 0) return(NULL)
  filtered.dots <- lapply(dots, function(x) x[keep])    
  
  return(c(list("alignments" = keep_alns), filtered.dots))
}

