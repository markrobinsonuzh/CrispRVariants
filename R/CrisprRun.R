#'@title CrisprRun class
#'@description A container for a single sample of alignments narrowed 
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
#'@seealso \code{\link[crispRvariants]{CrisprSet}}
#'@author Helen Lindsay
#'@export CrisprRun
#'@exportClass CrisprRun
CrisprRun = setRefClass(
  Class = "CrisprRun",
  fields = c(alns = "GAlignments",
             name = "character",
             query_ranges = "CompressedIRangesList",
             ref_ranges = "CompressedIRangesList",
             genome_ranges = "CompressedIRangesList",
             cigar_ops = "CompressedCharacterList",   
             insertions = "data.frame",
             ins_key = "integer", 
             cigar_labels = "character",
             chimeras = "GAlignments")
)

CrisprRun$methods(
  initialize = function(bam, target, genome.ranges, rc = FALSE, name = NULL, 
                        chimeras = GenomicAlignments::GAlignments(), verbose = TRUE){
    #Attributes:
    # cigar_labels are labels for variant combinations, e.g. used in plotting
    
    name <<- ifelse(is.null(name), NA, name)
    if (verbose == TRUE) cat(sprintf("\nInitialising CrisprRun %s\n", .self$name))
    
    alns <<- bam
    chimeras <<- chimeras 
    if (length(bam) == 0) { return() } 
    
    genome_ranges <<- genome.ranges
    ref_ranges <<- cigarRangesAlongReferenceSpace(cigar(bam))
    query_ranges <<- cigarRangesAlongQuerySpace(cigar(bam))
    
    cigar_ops <<- CharacterList(explodeCigarOps(cigar(bam)))
    .self$getInsertionSeqs()
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
    .self$field("query_ranges", .self$query_ranges[-idxs])
    .self$field("genome_ranges", .self$genome_ranges[-idxs])
    .self$field("ref_ranges", .self$ref_ranges[-idxs])
    .self$field("cigar_ops", .self$cigar_ops[-idxs])
  },

  countDeletions = function(count_multi_del = FALSE, count_del_w_ins = FALSE){
  '
Description:
Counts the number of reads containing a deletion

Input parameters:
  count_multi_del:   If TRUE, returns the exact number of deletions,
                     i.e., if one read contains 2 deletions, it contributes 2 to the
                     total count (default: FALSE)
  count_del_w_ins:   If TRUE, counts deletions regardless of whether reads also
                     contain insertions.  If FALSE, counts reads that contain 
                     deletions but not insertions (default: FALSE)
  ' 
    if (count_del_w_ins){
      if (count_multi_del)  return(sum(.self$cigar_ops %in% c("D", "N")))
      return(sum(any(.self$cigar_ops %in% c("D", "N"))))
    }
    if (count_multi_del){
      return(sum(.self$cigar_ops %in% c("D", "N")[! any.self$cigar_ops == "I"]))
    }
    return(sum(any(.self$cigar_ops %in% c("D", "N")) & ! any(.self$cigar_ops == "I")))
  },

  countInsertions = function(count_ins_w_del = FALSE, count_multi_ins = FALSE){
  '
Description:
Counts the number of reads containing an insertion

Input parameters:
  count_multi_ins:   If TRUE, returns the exact number of insertions,
                     i.e., if one read contains 2 insertions, it contributes 2 to the
                     total count (default: FALSE)
  count_ins_w_del:   If TRUE, counts insertions regardless of whether reads also
                     contain deletions  If FALSE, counts reads that contain 
                     insertions but not deletions (default: FALSE)
  '
    if (count_ins_w_del){
      if (count_multi_ins)  return(sum(.self$cigar_ops == "I"))
      return(sum(any(.self$cigar_ops == "I")))
    }
  
    if (count_multi_ins){
      return(sum(.self$cigar_ops == "I" & ! any(.self$cigar_ops %in% c("D","N"))))
    }
    return(sum(any(.self$cigar_ops == "I" & ! any(.self$cigar_ops %in% c("D","N")))))
  },
  
  countIndels = function(){
'
Description:
    Prints the number of target reads that include at least one 
    insertion or deletion, by counting cigar operations "I" (insertion),
    "D" (deletion) and "N" (junction operation used by some aligners)
'
    return(sum(any(.self$cigar_ops %in% c("I", "D", "N"))))
  },
  
  indelPercent = function(){
'
Description:
    Prints the percentage of target reads that include at least one 
    insertion or deletion
'
    return((.self$countIndels() / length(.self$cigar_ops))*100)
  },
  
  getInsertionSeqs = function(){ 
    # Note that the start of a ref_ranges insertion is its genomic end (rightmost base)
    
    ins <- .self$cigar_ops == "I"  
    idxs <- rep(1:length(.self$cigar_ops), sum(ins))
    tseqs <- as.character(mcols(.self$alns)$seq)[idxs]    
    
    if (length(tseqs) == 0) {
      insertions <<- data.frame()
      ins_key <<- integer()
      return()
    }    
    
    qranges <- unlist(.self$query_ranges[ins]) 
    ins_seqs <- as.character(subseq(tseqs, start(qranges), end(qranges)))
    ins_starts <- start(unlist(.self$ref_ranges[ins]))
    genomic_starts <- unlist(start(.self$genome_ranges[ins])) -1 # -1 for leftmost base
    
    df <- data.frame(start = ins_starts, seq = ins_seqs, genomic_start = genomic_starts)
    df$seq <- as.character(df$seq)
    insertions <<- aggregate(rep(1, nrow(df)), by = as.list(df), FUN = table)
    colnames(insertions) <<- c("start", "seq", "genomic_start", "count")
    # Store a key to match the sequences to their insertion
    ins_key <<- match(interaction(df), interaction(insertions[,c(1:3)]))
    names(ins_key) <<- idxs
  },
  
  getInsertionsTable = function(){
    ins_starts <- .self$insertions$genomic_start[.self$ins_key]
    
    if (length(ins_starts) == 0){
      ins <- data.frame()
    } else {
      ins <- data.frame(read = as.integer(names(.self$ins_key)), chromosome = chrom, 
                        start = ins_starts, end = ins_starts, 
                        inserted = .self$insertions$seq[.self$ins_key], stringsAsFactors = FALSE)
    }
    return(ins)
  },
  
  .checkNonempty = function(){
    if (length(.self$alns) == 0){
      message("No on target alignments")
      return(FALSE)
    }
    return(TRUE)
  },
  
  getVariants = function(ref_genome, chrom = NULL, ensembl = FALSE, strand = "+"){
'
Description:
  Returns a data frame of unique variants and their coordinates
    
Input parameters:
    ref_genome:   a BSGenome obj 
    chrom:        chromosomes to consider 
    ensembl:      should variants be returned in Ensembl default format, 
                  e.g. for use with the Ensembl Variant Effect Predictor 
                  (default: FALSE)
    strand:       The strand of the target alignments.  If "-", the 
                  reference sequence is reverse complemented and variants are 
                  returned w.r.t the negative strand (default: "+")
Result:
    Value "read" is the index of the alignment with the variant in .self$alns 
    '    
    if (! .self$.checkNonempty()){
      return(list())
    }
    
    # Get the reference sequence
    if (is.null(chrom)){
      chrom <- as.character(seqnames(.self$alns)[1])
    }
    if (! grepl("chr", chrom)){       
      chrom <- paste0("chr", chrom)
    }
    
    get_start <- min(start(.self$alns))
    get_end <- max(end(.self$alns))
    ref_seq <- getSeq(ref_genome, chrom, get_start, get_end)
    ref_offset <- start(.self$alns) - get_start
    
    if (strand == "-"){  # alignments and ranges have already been reverse complemented
      ref_seq <- reverseComplement(ref_seq)
      ref_offset <- get_end - end(.self$alns) 
    }
    
    # Get deletions:
    dels <- crun$cigar_ops %in% c("D", "N")
    del_idxs <- rep(seq_along(dels), sum(dels))
    
    if (length(del_idxs) == 0){
      dels <- data.frame()
    }else{
      del_ranges <- unlist(.self$ref_ranges[dels]) 
      del_seqs <- as.character(Views(ref_seq, shift(del_ranges, ref_offset[del_idxs])))
      del_gen <- unlist(.self$genome_ranges[dels])
      dels <- data.frame(read = del_idxs, chromosome = chrom, start = start(del_gen), 
                         end = end(del_gen), deleted = del_seqs, stringsAsFactors = FALSE)
    }
    
    # Get insertions: already computed at initialisation
    ins_starts <- .self$insertions$genomic_start[.self$ins_key]
    
    if (length(ins_starts) == 0){
      ins <- data.frame()
    } else {
      ins <- data.frame(read = as.integer(names(.self$ins_key)), chromosome = chrom, 
                        start = ins_starts, end = ins_starts, 
                        inserted = .self$insertions$seq[.self$ins_key], stringsAsFactors = FALSE)
    }
    
    # Get mismatches - check for mismatches and the read ranges with "M" (match/mismatch)
    matches <- lapply(.self$cigar_ops, function(x) which(x == "M"))
    match_ranges <- .self$query_ranges[matches]
    key <- rep(1:length(match_ranges), lapply(match_ranges, length))
    match_ranges <- unlist(match_ranges)
    match_seqs <- subseq(mcols(.self$alns)$seq[key], start(match_ranges), end(match_ranges))
    
    ref_ranges_ <- unlist(shift(.self$ref_ranges[matches], ref_offset))
    ref_seqs <- as(Views(ref_seq, ref_ranges_), "DNAStringSet")    
    idxs <- which(ref_seqs != match_seqs)  
    
    if (strand == "-"){
      # Make sequences read from genomic left to right for ease of getting genomic coords  
      match_seqs <- DNAStringSet(lapply(match_seqs, rev))
      ref_seqs <- DNAStringSet(lapply(ref_seqs, rev))
    }
    
    gen_starts <- start(unlist(.self$genome_ranges[matches]))
    
    mm <- do.call(rbind, lapply(idxs, function(i){
      r1 <- as.matrix(ref_seqs[i])
      m1 <- as.matrix(match_seqs[i])
      r1[m1 == "N"] <- "N"
      neq <- m1 != r1 
      if (length(which(neq)) == 0){ # Mismatches are only due to "N"
        return(data.frame())
      }
      gen_starts_i <- gen_starts[i] + (which(neq) -1 )
      df <- data.frame(read = key[i], chromosome = chrom, 
                       start =gen_starts_i, ref = r1[neq], query = m1[neq])
      
    }))
    
    vars <- list(insertions = ins, deletions = dels, mismatches = mm)  
    
    if (ensembl == TRUE){
      return(getVarsEnsemblFormat(vars, strand))
    } else{
      return(vars)
    }
  },
  
  getVarsEnsemblFormat = function(vars, strand = "*"){
    # Ensembl requires end = start - 1 for insertions
    # End points are included
    
    if (! .self$.checkNonempty()){
      return(list())
    }
    
    vins <- vars$insertions
    ins <- data.frame(chrom = vins$chromosome, start = vins$start, end = vins$end -1 , 
                      allele = sprintf("-/%s", vins$inserted), id = vins$read)
    
    vdels <- vars$deletions
    dels <- data.frame(chrom = vdels$chromosome, start = vdels$start, end = vdels$end,
                       allele = sprintf("%s/-", vdels$deleted), id = vdels$read)
    
    vmm <- vars$mismatches
    mm <- data.frame(chrom = vmm$chromosome, start = vmm$start, end = vmm$start,
                     allele = sprintf("%s/%s", vmm$ref, vmm$query), id = vmm$read)
    
    vars <- rbind(ins,dels,mm)
    if (nrow(vars) > 0){
      vars$strand <- strand
      vars <- vars[order(vars$id),]
      vars$id <- names(.self$alns[vars$id])
      vars <- vars[,c("chrom","start","end","allele","strand","id")]
    }
    vars$chrom <- gsub("chr", "", vars$chrom)
    vars
  },
  
  getCigarLabels = function(short = TRUE, match_label = "no variant", 
                            genome_to_target = NULL, rc = FALSE, target_start = NULL,
                            target_end = NULL, ref = NULL, split_non_indel = TRUE, 
                            mismatch_label = "SNV", cut.site = 18, upstream = 8, 
                            downstream = 5){
    
    # Sets / returns cigar labels                        
    # If short = only consider insertions and deletions (except for wildtype?)
    
    # genome_to_target, if provided, is a vector with names being genomic start coords
    # and values being coords with respect to the target site  
    
    # Shorten? -> Renumber?  If not renumbered, record starting location when different from 
    # target start
    if (length(.self$alns) == 0) return(vector())
    
    if (! length(.self$cigar_labels) == 0){
      return(.self$cigar_labels)
    }
    
    # for default cigar (short = FALSE) and short cigar without renumbering:  
    # indicate the starting location with respect to the target start if it is not the target
    start_offset <- rep("", length(.self$alns))
    
    if (short == FALSE | is.null(genome_to_target)){
      if (rc == TRUE){
        # Note: genome_ranges for rc are written rightmost to leftmost
        
        start_offset <- target_end - unlist(lapply(.self$genome_ranges, function(x) max(end(x))),
                                            use.names = FALSE)
        start_offset_new <- target_end - max(end(.self$genome_ranges))
        
      } else {
        start_offset <- unlist(lapply(.self$genome_ranges, function(x) min(start(x))), 
                               use.names = FALSE) - target_start
        start_offset_new <-  min(start(.self$genome_ranges)) - target_start 
      }
      is.nonzero <- start_offset != 0
      start_offset[is.nonzero] <- sprintf("(%s)", start_offset[is.nonzero])
      start_offset[! is.nonzero] <- ""       
      
      if ( short == FALSE ){
        cigar_labels <<- paste0(start_offset, cigar(.self$alns))
        return(.self$cigar_labels)
      }
    }
    
    #idxs <- .self$cigar_ops != "M"
    #ops <- .self$cigar_ops[idxs]
    #rranges <- .self$ref_ranges[idxs]
    #qranges <- .self$query_ranges[idxs]
    #cig_idxs <- rep(1:length(idxs), lapply(idxs, length))
    
    idxs <- lapply(.self$cigar_ops, function(x) which(x != "M"))
    ops <- lapply(seq_along(idxs), function(i) .self$cigar_ops[[i]][idxs[[i]]])
    rranges <- .self$ref_ranges[idxs]
    qranges <- .self$query_ranges[idxs]
    cig_idxs <- rep(1:length(idxs), lapply(idxs, length))   
    
    
    if (is.null(genome_to_target)){
      start <- unlist(start(rranges))
      end <- unlist(end(rranges))
    } else {
      gen_ranges <- .self$genome_ranges[idxs]     
      if (rc == FALSE){
        start <- genome_to_target[as.character(unlist(start(gen_ranges)))]
        end <- genome_to_target[as.character(unlist(end(gen_ranges)))]
        
      }else {
        start <- genome_to_target[as.character(unlist(end(gen_ranges)))]
        end <- genome_to_target[as.character(unlist(start(gen_ranges)))]
      }
    }
    
    info <- data.frame(op = unlist(ops), qwidth = unlist(width(qranges)), 
                       start = start, end = end, rwidth = unlist(width(rranges)))
    
    getShortOp <- function(mutn){
      if (mutn["op"] %in% c("N", "D")){ 
        sprintf('%s:%sD', as.numeric(mutn["start"]), as.numeric(mutn["rwidth"]))
      }else if (mutn["op"] == "I"){
        sprintf('%s:%sI', as.numeric(mutn["start"]), as.numeric(mutn["qwidth"]))
      }
      
    }
    result <- rep(match_label, length(idxs))
    if (nrow(info) > 0) {
      short_ops <- apply(info, 1, getShortOp)
      short_ops <- split(short_ops, cig_idxs)
      short_ops <- sapply(short_ops, function(x) do.call(paste, c(as.list(x),sep = ",")))
      result[as.numeric(names(short_ops))] <- short_ops
    }
    if (split_non_indel){
      result <- .splitNonIndel(ref, result, match_label, mismatch_label, cut.site, 
                               upstream, downstream)
    }
    cigar_labels <<- result
    result
  },
  
  .splitNonIndel = function(ref, cig_labels, match_label = "no variant", 
                          mismatch_label = "SNV", cut_site = 17, upstream = 8, 
                          downstream = min(5, width(ref) - cut_site)){
  
  # Only consider mismatches up to (upstream) to the left of the cut and
  # (downstream) to the right of the cut
  # The cut site is between cut_site and cut_site + 1
  
  # TO DO - AMBIGUITY CHARACTERS SHOULD NOT BE SNVS!  
    
  no_var <- which(cig_labels == match_label & mcols(.self$alns)$seq != ref)
  if (length(no_var) == 0) return(cig_labels)
  
  if ((cut_site-upstream + 1) < 0 | (cut_site + downstream) > length(ref)){
    stop("Specified range for detecting SNVs is greater than target range")    
  }
  
  snv_range <- c((cut_site-upstream + 1):(cut_site + downstream))
  no_var_seqs <- as.matrix(mcols(.self$alns)$seq[no_var])
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