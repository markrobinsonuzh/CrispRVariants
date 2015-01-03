#'@title CrisprSet class
#'@description A container for holding a set of narrowed alignments, 
#'each corresponding to the same target region.  Individual samples are 
#'represented as CrisprRun objects.  
#'@param crispr.runs A list of CrisprRun objects, typically representing individual samples
#'within an experiment
#'@param reference The reference sequence, must be the same length as the target region
#'@param target The target location (GRanges).  Variants will be counted over this region.
#'Need not correspond to the guide sequence.  
#'@param rc Should the alignments be reverse complemented, 
#'i.e. displayed w.r.t the reverse strand? (default: FALSE)
#'@param short_cigars If TRUE, variants labels are created from the location of their
#'insertions and deletions.  For variants with no insertions or deletions, the locations 
#'of any single base mismatches are displayed (default: TRUE).
#'@param names A list of names for each of the samples, e.g. for displaying in plots.
#'If not supplied, the names of the crispr.runs are used, which default to the filenames 
#'of the bam files if available (Default: NULL)
#'@param renumbered Should the variants be renumbered using target.loc as the zero point? 
#'If TRUE, variants are described by the location of their 5'-most base with respect to the 
#'target.loc.  A 3bp deletion starting 5bp 5' of the cut site would be labelled
#'(using short_cigars) as -5:3D (Default: TRUE)
#'@param target.loc The location of the Cas9 cut site with respect to the supplied target.
#'(Or some other central location).  Can be displayed on plots and used as the zero point 
#'for renumbering variants. For a target region with the PAM location from bases 21-23, 
#'the target.loc is base 18 (default: NA)
#'@param match.label Label for sequences with no variants (default: "no variant")
#'@param verbose If true, prints information about initialisation progress (default: TRUE)
#'@param ...
#'@field crispr_runs A list of CrisprRun objects, typically corresponding to samples 
#'of an experiment.  
#'@field ref The reference sequence for the target region, as a DNAString object 
#'@field cigar_freqs A matrix of counts for each variant
#'@field target The target location, as a GRanges object
#'@author Helen Lindsay
#'@seealso \code{\link[crispRvariants]{CrisprRun}}
#'@export CrisprSet
#'@exportClass CrisprSet 
CrisprSet = setRefClass(
  Class = "CrisprSet",
  fields = c(crispr_runs = "list", 
             ref = "DNAString",
             insertion_sites = "data.frame",
             cigar_freqs = "matrix",
             target = "GRanges",
             genome_to_target = "integer",
             pars = "list")
)

CrisprSet$methods(
  initialize = function(crispr.runs, reference, target, rc = FALSE, short_cigars = TRUE, 
                        names = NULL, renumbered = TRUE, target.loc = NA, 
                        match.label = "no variant", verbose = TRUE, ...){
    
    print(sprintf("Initialising CrisprSet with %s samples", length(crispr.runs)))
    
    reference <- as(reference, "DNAString")
    if (width(target) != length(reference)){
      stop("The target and the reference sequence must be the same width")
    }
    if (renumbered == TRUE & is.na(target.loc)){
      stop(paste0("Must specify target.loc for renumbering variant locations.\n",
                  "The target.loc is the zero point with respect to the reference string.\n",
                  "This is typically 18 for a 23 bp Crispr-Cas9 guide sequence"))
    }
    
    # TO DO - DECIDE WHETHER TO KEEP EMPTY RUNS
    target <<- target   
    ref <<- reference 
    pars <<- list("match_label" = match.label, "target.loc" = target.loc, 
                  "mismatch_label" = "SNV", "renumbered" = renumbered)
    
    #pars <<- modifyList(pars, ...)   
    
    crispr_runs <<- crispr.runs
    
    if (is.null(names)) {
      names(.self$crispr_runs) <- sapply(.self$crispr_runs, function(x) x$name)
    }else {
      names(.self$crispr_runs) <- names
    }
    nonempty_runs <- sapply(.self$crispr_runs, function(x) {
      ! length(x$alns) == 0})
    
    .self$crispr_runs <<- .self$crispr_runs[nonempty_runs]
    if (length(.self$crispr_runs) == 0) stop("no on target runs")
    
    if (verbose == TRUE) cat("Renaming cigar strings and counting indels\n")
    cig_by_run <- .self$.setCigarLabels(renumbered = renumbered, target.loc = target.loc,
                                        target_start = start(target), target_end = end(target), 
                                        rc = rc, match_label = match.label, ref = ref)
    .self$.countCigars(cig_by_run)
  },
  
  show = function(){
    print(sprintf(paste0("CrisprSet object containing %s CrisprRun samples\n", 
                         "Target location:\n"), length(.self$crispr_runs)))
    print(.self$target)
    print("Most frequent variants:")
    print(.self$.getFilteredCigarTable(top_n = 6))
  },
  
  .setCigarLabels = function(renumbered = FALSE, target.loc = NA, target_start = NA,
                             target_end = NA, rc = FALSE, match_label = "no variant",
                             ref = NULL){
    g_to_t <- NULL
    
    if (renumbered == TRUE){
      if (any(is.na(c(target_start, target_end, rc)))){
        stop("Must specify target.loc (cut site), target_start, target_end and rc
             for renumbering")
      }
      g_to_t <- genomeToTargetLocs(target.loc, target_start, target_end, rc)
      }
    
    cut.site <- ifelse(is.na(target.loc), 18, target.loc)
    
    cig_by_run <- lapply(.self$crispr_runs,
                         function(crun) crun$getCigarLabels(match_label = match_label, rc = rc,
                                                            genome_to_target = g_to_t, ref = ref, 
                                                            cut.site = cut.site))     
    return(cig_by_run)
  }, 
  
  .countCigars = function(cig_by_run = NULL, nonvar_first = TRUE){
    # Note that this function does not consider starts, two alignments starting at
    # different locations but sharing a cigar string are considered equal
    
    if (is.null(cig_by_run)){
      cig_by_run <- lapply(.self$crispr_runs, function(crun) crun$getCigarLabels())
    }
    
    unique_cigars <- unique(unlist(cig_by_run))  
    
    m <- matrix(unlist(lapply(cig_by_run, function(x) table(x)[unique_cigars])), 
                nrow = length(unique_cigars), dimnames = list(unique_cigars, names(.self$crispr_runs)))
    
    m[is.na(m)] <- 0
    m <- m[order(rowSums(m), decreasing = TRUE),, drop = FALSE]
    
    if (nonvar_first){    
      is_ref <- grep(.self$pars["match_label"], rownames(m))
      is_snv <- grep(.self$pars["mismatch_label"], rownames(m))
      new_order <- setdiff(1:nrow(m), c(is_ref,is_snv))
      m <- m[c(is_ref, is_snv, setdiff(1:nrow(m), c(is_ref,is_snv))),,drop = FALSE]
    }
    
    cigar_freqs <<- m
  },
  
  .getFilteredCigarTable = function(top_n = nrow(.self$cigar_freqs), freq_cutoff = 0){
    rs <- rowSums(.self$cigar_freqs)
    minfreq <- rs >= freq_cutoff
    topn <- rank(-rs) <= top_n
    cig_freqs <- .self$cigar_freqs[minfreq & topn ,, drop = FALSE] 
    return(cig_freqs) 
  },
  
  countVariantAlleles = function(counts_t = NULL){
    # Returns counts of variant alleles
    # SNV alleles are not considered variants here
    if (is.null(counts_t)) counts_t <- .self$cigar_freqs
    counts_t <- counts_t[!rownames(counts_t) == .self$pars["match_label"],,drop = FALSE]
    alleles <- colSums(counts_t != 0)    
    return(data.frame(Allele = alleles, Sample = names(alleles)))
  },
  
  heatmapCigarFreqs = function(as_percent = FALSE, x_size = 16, y_size = 16, 
                               x_axis_title = NULL, x_angle = 90, annotate_counts = TRUE, 
                               freq_cutoff = 0, top_n = nrow(.self$cigar_freqs), ...){
    
    cig_freqs <- .getFilteredCigarTable(top_n, freq_cutoff)
    p <- cigarFrequencyHeatmap(cig_freqs, as_percent, x_size, y_size, x_axis_title,
                               x_angle, annotate_counts, ...)
    return(p)
  },
  
  plotVariants = function(freq_cutoff = 0, top_n = nrow(.self$cigar_freqs), 
                          short_cigars = FALSE, renumbered = .self$pars["renumbered"], 
                          show_genomic = FALSE, ...){
    
    # var_freq_cutoff = i (integer) only plot variants that occur >= i times
    # top_n = total number of variants to plot
    # ... arguments for plotAlignments
    # note that if there are ties, top_n only includes ties with 
    # all members ranking <= top_n 
    
    cig_freqs <- .self$.getFilteredCigarTable(top_n, freq_cutoff)
    
    alns <- .self$makePairwiseAlns(cig_freqs)
    if (!("cigar" %in% colnames(.self$insertion_sites))){
      .self$getInsertions() 
    }
    
    # How should the x-axis be numbered? 
    # Baseline should be numbers, w optional genomic locations
    if (renumbered == TRUE){
      genomic_coords <- c(start(.self$target):end(.self$target))
      target_coords <- .self$genome_to_target[as.character(genomic_coords)]
      if (as.character(strand(.self$target)) == "-"){
        target_coords <- rev(target_coords)
      }
      xbreaks = which(target_coords %% 5 == 0 | abs(target_coords) == 1)
      target_coords <- target_coords[xbreaks]
      
      p <- plotAlignments(.self$ref, alns, .self$insertion_sites, 
                          xtick_labs = target_coords, xtick_breaks = xbreaks, ...)
    } else {
      p <- plotAlignments(.self$ref, alns, .self$insertion_sites, ...)    
    }
    return(p)
  },
  
  getInsertions = function(with_cigars = TRUE){
    if (with_cigars == FALSE){
      all_ins <- do.call(rbind, lapply(.self$crispr_runs, function(x) x$insertions))
    } else {
      all_ins <- do.call(rbind, lapply(.self$crispr_runs, function(x) {
        ik <- x$ins_key
        v <- data.frame(ik, x$getCigarLabels()[as.integer(names(ik))])
        v <- v[!duplicated(v),]
        v <- v[order(v$ik),]
        cbind(x$insertions[v[,1],], cigar = v[,2])
      }))
    }
    if (nrow(all_ins) == 0) {
      insertion_sites <<- all_ins
      return()
    }
    insertion_sites <<- all_ins[order(all_ins$start, all_ins$seq),, drop = FALSE]
  },
  
  makePairwiseAlns = function(cig_freqs = .self$cigar_freqs, ...){
    # Get alignments by cigar string, make the alignment for the consensus
    # The short cigars (not renumbered) do not have enough information, 
    # use the full cigars for sorting
    
    cigs <- unlist(lapply(.self$crispr_runs, function(x) cigar(x$alns)), use.names = FALSE)
    cig_labels <- unlist(lapply(.self$crispr_runs, function(x) x$getCigarLabels()), use.names = FALSE)
    
    names(cigs) <- cig_labels # calling by name with duplicates returns the first match
    
    splits <- split(seq_along(cig_labels), cig_labels)
    splits <- splits[match(rownames(cig_freqs), names(splits))]
    
    splits_labels <- names(splits)
    names(splits) <- cigs[names(splits)]
    
    x <- lapply(.self$crispr_runs, function(x) x$alns)
    all_alns <- do.call(c, unlist(x, use.names = FALSE))
    
    seqs <- c()
    starts <- c()
    
    # SOMEWHERE HERE - CONSENSUS SEQ ONLY TAKING ONE SAMPLE?
    
    for (i in seq_along(splits)){
      idxs <- splits[[i]]
      seqs[i] <- consensusString(mcols(all_alns[idxs])$seq)
      start <- unique(start(all_alns[idxs]))
      if (length(start) > 1)
        stop("Sequences with the same cigar string have different starting locations.
             This case is not implemented yet.")  
      starts[i] <- start[1]
    }  
    
    alns <- mapply(seqsToAln, names(splits), seqs, aln_start = starts, 
                   target_start = start(.self$target), target_end = end(.self$target), ...)
    
    names(alns) <- splits_labels
    alns
  },
  
  plotVariantOverview = function(){
    
    # heatmap must take same filtering args, OR be able to match the alignment names
    
  },
  
  genomeToTargetLocs = function(target.loc, target_start, target_end, rc = FALSE){
    # target.loc should be relative to the start of the target sequence, even if the 
    # target is on the negative strand
    # target.loc is the left side of the cut site (Will be numbered -1)
    # target_start and target_end are genomic coordinates, with target_start < target_end
    # rc: is the target on the negative strand wrt the reference?
    # returns a vector of genomic locations and target locations
    
    # Example:  target.loc = 5
    # Before: 1  2  3  4  5  6  7  8 
    # After: -5 -4 -3 -2 -1  1  2  3
    # Left =  original - target.loc - 1
    # Right = original - target.loc
    
    gs <- min(unlist(lapply(.self$crispr_runs, function(x) start(unlist(x$genome_ranges)))))
    
    all_gen_ranges <- lapply(.self$crispr_runs, function(x) x$genome_ranges)
    
    # Nope - here gives a vector
    #gs_test <- min(start(do.call(c, unlist(all_gen_ranges, use.names = FALSE))))
    #cat(sprintf("testing gs: %s = %s: %s", gs_test, gs, gs_test == unlist(gs, use.names=FALSE)))
    
    ge <- max(unlist(lapply(.self$crispr_runs, function(x) end(unlist(x$genome_ranges)))))
    
    if (rc == TRUE){
      tg <- target_end - (target.loc - 1)
      new_numbering <- rev(c(seq(-1*(ge - (tg -1)),-1), c(1:(tg - gs))))
      names(new_numbering) <- c(gs:ge)
      
    } else {
      tg <- target_start + target.loc - 1
      new_numbering <- c(seq(-1*(tg - (gs-1)),-1), c(1:(ge - tg)))
      names(new_numbering) <- c(gs:ge)
    }
    genome_to_target <<- new_numbering
    new_numbering
  }
)