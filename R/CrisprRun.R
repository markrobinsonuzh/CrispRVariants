#'@title CrisprRun class
#'@description A ReferenceClass container for a single sample of alignments narrowed 
#'to a target region
#'@param bam a GAlignments object containing (narrowed) alignments to the target region.
#' Filtering of the bam should generally be done before initialising a CrisprRun object
#'@param target The target location, a GRanges object
#'@param genome.ranges A GRangesList of genomic coordinates for the cigar operations.
#' If bam is a standard GAlignments object, this is equivalent to 
#' cigarRangesAlongReferenceSpace + start(bam)
#'@param rc (reverse complement)  Should the alignments be reverse complemented,
#'i.e. displayed with respect to the negative strand?  (Default: FALSE)
#'@param name A name for this set of reads, used in plots if present (Default: NULL)
#'@param chimeras Off-target chimeric alignments not in bam.  (Default: empty)
#'@param verbose Print information about initialisation progress (Default: TRUE) 
#'@field alns A GAlignments object containing the narrowed reads.  Note that if the alignments
#'are represented with respect to the reverse strand, the "start" remains with repect to the
#'forward strand, whilst the cigar and the sequence are reverse complemented.  
#'@field name The name of the sample
#'@field cigar_labels A vector of labels for the reads, based on the cigar strings,
#'optionally renumbered with respect to a new zero point (e.g. the cut site) and 
#'shortened to only insertion and deletion locations.  
#'Set at initialisation of a CrisprSet object, but not at 
#'initialisation of a CrisprRun object.
#'@field chimeras Chimeric, off-target alignments corresponding to alignments in alns
#'@seealso \code{\link[CrispRVariants]{CrisprSet}}
#'@author Helen Lindsay
#'@export CrisprRun
#'@exportClass CrisprRun
CrisprRun = setRefClass(
  Class = "CrisprRun",
  fields = c(alns = "GAlignments",
             name = "character",
             cigar_ops = "CompressedCharacterList",   
             insertions = "data.frame",
             ins_key = "integer", 
             cigar_labels = "character",
             chimeras = "GAlignments",
             chimera_combs = "data.frame")
)

CrisprRun$methods(
  initialize = function(bam, target, genome.ranges, rc = FALSE, name = NULL, 
                        chimeras = GenomicAlignments::GAlignments(), verbose = TRUE){
    
    name <<- ifelse(is.null(name), NA, name)
    if (verbose == TRUE) cat(sprintf("\nInitialising CrisprRun %s\n", .self$name))
    
    alns <<- bam
    chimeras <<- chimeras 
    chimera_combs <<- .self$splitChimeras()
    
    if (length(bam) == 0) { return() } 
    
    ref_ranges <- cigarRangesAlongReferenceSpace(cigar(bam))
    cigar_ops <<- CharacterList(explodeCigarOps(cigar(bam)))
    
    # recalculate this in case of keep_unpaired = FALSE, not tested
    genome.ranges <- cigarRangesAlongReferenceSpace(cigar(bam), pos = start(bam))
    .self$getInsertionSeqs(ref_ranges = ref_ranges, genome_ranges = genome.ranges)
  },
  
  show = function(){
    print(c(class(.self), sprintf("CrisprRun object named %s, with %s on target alignments.", 
                                  .self$name, length(.self$alns)), .self$alns))
  },
  
  removeSeqs = function(idxs){
'
Description:
  Remove sequences from a CrisprRun object and from the internal CrisprRun 
  fields that store insertion locations for plotting.

Input parameters:
  idxs:     Indexes of reads to remove'
    
    # note insertions table is not updated
    ins_key_idxs <- which(names(.self$ins_key) %in% idxs)
    
    if (length(ins_key_idxs) > 0){
      .self$field("ins_key", .self$ins_key[-ins_key_idxs]) 
    }
    
    # Insertion key refers to indexs with in the alignments, these have shifted after 
    # filtering.  Note - would need to shift insertions idxs if insertions is also updated
    subtract <- rep(0, length(.self$alns))
    subtract[idxs] <- 1
    subtract <- cumsum(subtract)
    
    temp <- .self$ins_key 
    nm_as_num <- as.numeric(names(temp))
    names(temp) <- nm_as_num - subtract[nm_as_num]
    .self$field("ins_key", temp)  
    
    # Remove the extra sequences from the chimeras
    rm_by_nm <- names(alns)[idxs]
    ch_to_keep <- !(names(.self$chimeras) %in% rm_by_nm)
    .self$field("chimeras", .self$chimeras[ch_to_keep])

    .self$field("alns", .self$alns[-idxs])
    .self$field("cigar_labels", .self$cigar_labels[-idxs])
    .self$field("cigar_ops", .self$cigar_ops[-idxs])
  },

  getInsertionSeqs = function(ref_ranges, genome_ranges){ 
    # Note that the start of a ref_ranges insertion is its genomic end (rightmost base)    
    
    ins <- .self$cigar_ops == "I"  
    idxs <- rep(1:length(.self$cigar_ops), sum(ins))
    tseqs <- as.character(mcols(.self$alns)$seq)[idxs]    
    
    if (length(tseqs) == 0) {
      insertions <<- data.frame()
      ins_key <<- integer()
      return()
    }    
    
    query_ranges <- cigarRangesAlongQuerySpace(cigar(.self$alns))
    qranges <- unlist(query_ranges[ins]) 
    ins_seqs <- as.character(subseq(tseqs, start(qranges), end(qranges)))
    ins_starts <- start(unlist(ref_ranges[ins]))
    genomic_starts <- unlist(start(genome_ranges[ins])) -1 # -1 for leftmost base
    df <- data.frame(start = ins_starts, seq = ins_seqs, genomic_start = genomic_starts)
    df$seq <- as.character(df$seq)
    insertions <<- aggregate(rep(1, nrow(df)), by = as.list(df), FUN = table)
    colnames(insertions) <<- c("start", "seq", "genomic_start", "count")
    # Store a key to match the sequences to their insertion
    ins_key <<- match(interaction(df), interaction(insertions[,c(1:3)]))
    names(ins_key) <<- idxs
  },
  
  .checkNonempty = function(){
    if (length(.self$alns) == 0){
      message("No on target alignments")
      return(FALSE)
    }
    return(TRUE)
  },
  
  splitChimeras = function(){
    splits <- split(cigar(.self$chimeras), names(.self$chimeras))
    if (length(splits) == 0) return(data.frame())
    combination <- sapply(splits, paste, collapse=";")
    tt <- as.data.frame(table(combination))
    tt <- tt[order(tt$Freq, decreasing = TRUE),]
    return(tt)
  },  


  
  getCigarLabels = function(target.loc, genome_to_target, ref,
                             separate.snv = TRUE, match.label = "no variant",
                             mismatch.label = "SNV",
                             rc = FALSE, keep.ops = c("I","D","N"), upstream = 8, 
                            downstream = min(5, width(ref) - cut_site)){
    
    cigs <- cigar(.self$alns)
    wdths <- explodeCigarOpLengths(cigs)
    ops <- explodeCigarOps(cigs)
    temp <- CharacterList(relist(paste0(unlist(wdths), unlist(ops)), wdths))
    ops <- CharacterList(ops)  
    keep <- ops %in% keep.ops
  
    rranges <- cigarRangesAlongReferenceSpace(cigs, pos = start(.self$alns))[keep]
    # here get snvs, add snv ops
  
    if (rc == TRUE){
      glocs <- end(rranges)
    } else {
      glocs <- start(rranges)
    }  
    temp <- paste(genome_to_target[as.character(unlist(glocs))],
                  unlist(temp[keep]), sep = ":")
    temp <- as.list(relist(temp, IRanges::PartitioningByEnd(cumsum(sum(keep)))))
    complex <- sum(keep) > 1

    if (rc == TRUE){
      temp[complex] <- sapply(temp[complex], function(x) paste(rev(x), collapse = ","))
    } else {
      temp[complex] <- sapply(temp[complex], paste, collapse = ",")
    }
    temp[sum(keep) == 0] <- match.label
    renamed <- as.character(temp)
    renamed <- .self$.splitNonIndel(ref, renamed, rc, match_label = match.label, 
                                    mismatch_label = mismatch.label, cut_site = target.loc,
                                    upstream = upstream, downstream = downstream)
    
    .self$field("cigar_labels", renamed)
    return(renamed)
  },

  .splitNonIndel = function(ref, cig_labels, rc, match_label = "no variant", 
                          mismatch_label = "SNV", cut_site = 17, 
                          upstream = 8, downstream = 5){
  
  # Only consider mismatches up to (upstream) to the left of the cut and
  # (downstream) to the right of the cut
  # The cut site is between cut_site and cut_site + 1
  upstream = min(cut_site, upstream)   
  downstream = min(downstream, nchar(ref) - cut_site)
  no_var <- which(cig_labels == match_label & mcols(.self$alns)$seq != ref)
  if (length(no_var) == 0) return(cig_labels)
  
  if ((cut_site-upstream + 1) < 0 | (cut_site + downstream) > length(ref)){
    stop("Specified range for detecting SNVs is greater than target range")    
  }
  
  snv_range <- c((cut_site-upstream + 1):(cut_site + downstream))
  sqs <- mcols(.self$alns)$seq[no_var]
  if (rc == TRUE) sqs <- reverseComplement(sqs)

  no_var_seqs <- as.matrix(sqs)
  no_var_seqs <- no_var_seqs[,snv_range, drop = FALSE]    
  
  rr <- strsplit(as.character(ref[snv_range]), "")[[1]]
  result <- apply(no_var_seqs, 1, function(x){
    snvs <- which((x != rr & x != "N")) - upstream - 1
    snvs[snvs >= 0] <- snvs[snvs >= 0] + 1 
    sprintf("%s:%s", mismatch_label, paste(snvs, collapse = ","))
  })
  result[result == sprintf("%s:", mismatch_label)] <- match_label
  cig_labels[no_var] <- result 
  return(cig_labels)       
  }
  
)