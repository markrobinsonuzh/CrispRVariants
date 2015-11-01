#'@title CrisprSet class
#'@description A ReferenceClass container for holding a set of narrowed alignments,
#'each corresponding to the same target region.  Individual samples are
#'represented as CrisprRun objects.  CrisprRun objects with no on-target
#'reads are excluded.
#'CrisprSet objects are constructed with \code{\link{readsToTarget}} or
#'\code{\link{readsToTargets}}.  For most use cases, a CrisprSet object should not
#'be accessed directly.
#'@param crispr.runs A list of CrisprRun objects, typically representing individual samples
#'within an experiment
#'@param reference The reference sequence, must be the same length as the target region
#'@param target The target location (GRanges).  Variants will be counted over this region.
#'Need not correspond to the guide sequence.
#'@param rc Should the alignments be reverse complemented,
#'i.e. displayed w.r.t the reverse strand? (default: FALSE)
#'@param short.cigars If TRUE, variants labels are created from the location of their
#'insertions and deletions.  For variants with no insertions or deletions, the locations
#'of any single base mismatches are displayed (default: TRUE).
#'@param names A list of names for each of the samples, e.g. for displaying in plots.
#'If not supplied, the names of the crispr.runs are used, which default to the filenames
#'of the bam files if available (Default: NULL)
#'@param renumbered Should the variants be renumbered using target.loc as the zero point?
#'If TRUE, variants are described by the location of their 5'-most base with respect to the
#'target.loc.  A 3bp deletion starting 5bp 5' of the cut site would be labelled
#'(using short.cigars) as -5:3D (Default: TRUE)
#'@param target.loc The location of the Cas9 cut site with respect to the supplied target.
#'(Or some other central location).  Can be displayed on plots and used as the zero point
#'for renumbering variants. For a target region with the PAM location from bases 21-23,
#'the target.loc is base 17 (default: NA)
#'@param match.label Label for sequences with no variants (default: "no variant")
#'@param mismatch.label Label for sequences with only single nucleotide variants
#'  (default: "SNV")
#'@param split.snv Should single nucleotide variants (SNVs) be shown for
#' reads without an insertion or deletion? (default: TRUE)
#'@param upstream.snv  If split.snv = TRUE, how many bases upstream of the target.loc
#' should SNVs be shown?  (default: 8)
#'@param downstream.snv If split.snv = TRUE, how many bases downstream of the target.loc
#' should SNVs be shown? (default: 5)
#'@param verbose If true, prints information about initialisation progress (default: TRUE)
#'@field crispr_runs A list of CrisprRun objects, typically corresponding to samples
#'of an experiment.
#'@field ref The reference sequence for the target region, as a Biostrings::DNAString object
#'@field cigar_freqs A matrix of counts for each variant
#'@field target The target location, as a GRanges object
#'@author Helen Lindsay
#'@seealso \code{\link{readsToTarget}} and \code{\link{readsToTargets}}
#'for initialising a CrisprSet, \code{\link[CrispRVariants]{CrisprRun}}
#'@examples
#'# Load the metadata table
#'md_fname <- system.file("extdata", "gol_F1_metadata_small.txt", package = "CrispRVariants")
#'md <- read.table(md_fname, sep = "\t", stringsAsFactors = FALSE)
#'
#'# Get bam filenames and their full paths
#'bam_fnames <- sapply(md$bam.filename, function(fn){
#'  system.file("extdata", fn, package = "CrispRVariants")})
#'
#'reference <- Biostrings::DNAString("GGTCTCTCGCAGGATGTTGCTGG")
#'gd <- GRanges("18", IRanges(4647377, 4647399), strand = "+")
#'
#'crispr_set <- readsToTarget(bam_fnames, target = gd, reference = reference,
#'                            names = md$experiment.name, target.loc = 17)
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
  initialize = function(crispr.runs, reference, target, rc = FALSE, short.cigars = TRUE,
                        names = NULL, renumbered = TRUE, target.loc = NA,
                        match.label = "no variant", mismatch.label = "SNV",
                        split.snv = TRUE, upstream.snv = 8, downstream.snv = 5,
                        verbose = TRUE, ...){

    if (isTRUE(verbose)){
      message(sprintf("Initialising CrisprSet %s:%s-%s with %s samples",
                  as.character(seqnames(target)), start(target), end(target),
                  length(crispr.runs)))
    }
    if (class(reference) == "DNAStringSet" || class(reference) == "character"){
      if (length(reference) > 1){
        stop("A CrisprSet contains alignments to exactly one reference sequence")
      }
      reference <- as(reference[[1]], "DNAString")
    }

    if (width(target) != length(reference)){
      stop("The target and the reference sequence must be the same width")
    }
    if (isTRUE(renumbered) & is.na(target.loc)){
      stop("Must specify target.loc for renumbering variant locations.\n",
           "The target.loc is the zero point with respect to the reference string.\n",
           "This is typically 17 for a 23 bp Crispr-Cas9 guide sequence")
    }

    target <<- target
    ref <<- reference
    pars <<- list("match_label" = match.label, "target.loc" = target.loc,
                  "mismatch_label" = mismatch.label, "renumbered" = renumbered,
                  "all_chimeric" = FALSE)

    crispr_runs <<- crispr.runs

    if (is.null(names)) {
      names(.self$crispr_runs) <- sapply(.self$crispr_runs, function(x) x$name)
    }else {
      names(.self$crispr_runs) <- names
    }
    nonempty_runs <- sapply(.self$crispr_runs, function(x) {
      ! ( length(x$alns) == 0 & length(x$chimeras) == 0 )
      })

    .self$crispr_runs <<- .self$crispr_runs[nonempty_runs]
    if (length(.self$crispr_runs) == 0) stop("no on target reads in any sample")

    if (unique( sapply(.self$crispr_runs, function(x) length(x$alns)) == 0)){
      pars["all_chimeric"] <<- TRUE
      return()
    }

    if (isTRUE(verbose)) message("Renaming cigar strings\n")

    cig_by_run <- .self$.setCigarLabels(renumbered = renumbered, target.loc = target.loc,
                                        target_start = start(target), target_end = end(target),
                                        rc = rc, match_label = match.label,
                                        mismatch_label = mismatch.label, ref = ref,
                                        short = short.cigars, split.snv = split.snv)
    if (isTRUE(verbose)) message("Counting variant combinations\n")
    .self$.countCigars(cig_by_run)
    .self$.getInsertions()
  },

  show = function(){
    cat(sprintf(paste0("CrisprSet object containing %s CrisprRun samples\n",
                         "Target location:\n"), length(.self$crispr_runs)))
    print(.self$target)
    print("Most frequent variants:")
    print(.self$.getFilteredCigarTable(top.n = 6))
  },

  .setCigarLabels = function(renumbered = FALSE, target.loc = NA, target_start = NA,
                             target_end = NA, rc = FALSE, match_label = "no variant",
                             mismatch_label = "SNV", short = TRUE, split.snv = TRUE,
                             upstream.snv = 8, downstream.snv = 5, ref = NULL){
    g_to_t <- NULL

    if (isTRUE(renumbered)){
      if (any(is.na(c(target_start, target_end, rc)))){
        stop("Must specify target.loc (cut site), target_start, target_end and rc
             for renumbering")
      }
      g_to_t <- .self$.genomeToTargetLocs(target.loc, target_start, target_end, rc)
      }
    cut.site <- ifelse(is.na(target.loc), 17, target.loc)

    # This section is slow
    cig_by_run <- lapply(.self$crispr_runs,
                             function(crun)  crun$getCigarLabels(
                              target.loc = cut.site, genome_to_target = g_to_t,
                              ref = .self$ref,
                              separate.snv = split.snv,
                              match.label = .self$pars$match_label,
                              mismatch.label = .self$pars$mismatch_label,
                              rc = rc, upstream = upstream.snv,
                              downstream = downstream.snv))

    cig_by_run
  },

  .countCigars = function(cig_by_run = NULL, nonvar_first = TRUE){
    # Note that this function does not consider starts, two alignments starting at
    # different locations but sharing a cigar string are considered equal

    if (is.null(cig_by_run)){
      cig_by_run <- lapply(.self$crispr_runs, function(crun) crun$cigar_label)
    }

    unique_cigars <- unique(unlist(cig_by_run))

    chimera_combs <- lapply(.self$crispr_runs, function(crun) crun$chimera_combs)

    m <- matrix(unlist(lapply(cig_by_run, function(x) table(x)[unique_cigars])),
                nrow = length(unique_cigars),
                dimnames = list(unique_cigars, names(.self$crispr_runs)))

    m[is.na(m)] <- 0
    m <- m[order(rowSums(m), decreasing = TRUE),, drop = FALSE]

    if (nonvar_first){
      is_ref <- grep(.self$pars["match_label"], rownames(m))
      is_snv <- grep(.self$pars["mismatch_label"], rownames(m))
      new_order <- setdiff(1:nrow(m), c(is_ref,is_snv))
      m <- m[c(is_ref, is_snv, setdiff(1:nrow(m), c(is_ref,is_snv))),,drop = FALSE]
    }

    .self$field("cigar_freqs", m)
  },

  filterUniqueLowQual = function(min_count = 2, max_n = 0, verbose = TRUE){
'
Description:
  Deletes reads containing rare variant combinations and more than
  a minimum number of ambiguity characters within the target region.
  These are assumed to be alignment errors.

Input parameters:
  min_count:    the number of times a variant combination must occur across
                all samples to keep (default: 2, i.e. a variant must occur
                at least twice in one or more samples to keep)
  max_n:        maximum number of ambiguity ("N") bases a read with a rare
                variant combination may contain.  (default: 0)
  verbose:      If TRUE, print the number of sequences removed (default: TRUE)
'
    # Find low frequency variant combinations, then find the corresponding samples
    low_freq <- .self$cigar_freqs[rowSums(.self$cigar_freqs) < min_count, , drop = FALSE]
    lf_cig_by_run <- apply(low_freq, 2, function(x) names(x)[x != 0])
    lns <- lapply(lf_cig_by_run, length)
    lf_cig_by_run <- lf_cig_by_run[lns > 0]

    # Find the corresponding reads, count the ambiguity characters
    rm_cset <- unlist(lapply(names(lf_cig_by_run), function(name){
      crun <- cset$crispr_runs[[name]]
      get_idxs <- match(lf_cig_by_run[[name]], crun$cigar_labels)
      sqs <- mcols(crun$alns[get_idxs])$seq
      to_remove <- get_idxs[as.numeric(Biostrings::letterFrequency(sqs, "N")) > max_n]
      cset_to_remove <- match(crun$cigar_labels[to_remove], rownames(cset$cigar_freqs))
      if (length(to_remove) > 0) crun$removeSeqs(to_remove)
      return(cset_to_remove)
    }))
    if (verbose){
      message(sprintf("Removing %s rare sequence(s) with ambiguities\n",
                      length(rm_cset)))
    }
    if ( length(rm_cset) > 0){
      .self$field("cigar_freqs", .self$cigar_freqs[-rm_cset,,drop = FALSE])
    }
  },

  .getFilteredCigarTable = function(top.n = nrow(.self$cigar_freqs),
                                    min.count = 1, min.freq = 0,
                                    include.chimeras = TRUE,
                                    include.nonindel = TRUE,
                                    type = c("counts", "proportions")){

    # Consider reordering by proportion instead of count at initialisation?

    result <- match.arg(type)

    # Add the chimeric alignments to the bottom
    if (isTRUE(include.chimeras)){
      ch_cnts <- sapply(.self$crispr_runs, function(crun) {
        length(unique(names(crun$chimeras)))
      })

      m <- .self$cigar_freqs
      if (length(m) > 0 & sum(ch_cnts) > 0){
        m_nms <- rownames(m)
        m <- rbind(m, ch_cnts)
        rownames(m) <- c(m_nms, "Other")
      } else if (length(m) == 0){
        m <- matrix(ch_cnts, nrow = 1, dimnames =list("Other", names(ch_cnts)))
      }
      # If all rows should be returned, add 1 to top.n
      if (top.n == nrow(.self$cigar_freqs)) top.n <- top.n + 1
    } else {
      m <- .self$cigar_freqs
    }

    # Filtering takes precedence over removing nonvariants
    # and selecting top.n

    # Default freq cutoff drops "Other" if there are no chimeras
    propns <- prop.table(m,2) * 100
    # At least one column has proportion greater than cutoff
    keep_freq <- rowSums(propns >= min.freq) > 0
    keep_count <- rowSums(m >= min.count) > 0
    propns <- propns[keep_freq & keep_count,, drop = FALSE]

    # Top variants are calculated by proportional contributions
    rs <- rowSums(propns)
    topn <- rank(-rs) <= top.n

    if (result == "proportions"){
      # Return proportions
      m <- propns[topn ,, drop = FALSE]
    } else {
      # Return counts
      m <- m[keep_freq & keep_count,, drop = FALSE][topn,, drop = FALSE]
    }

    if (include.nonindel == FALSE){
      nvr <- sprintf("%s|%s", .self$pars$match_label, .self$pars$mismatch_label)
      m <- m[!grepl(nvr, rownames(m)),,drop = FALSE]
    }

    m
  },

  .getUniqueIndelRanges = function(add_chr = TRUE, add_to_ins = TRUE){
    # Note this only gets the ranges, not the sequences, inserted sequences may differ
    # Returns a GRanges object of all insertions and deletions, with names = variant names
    # if "add_chr" == TRUE, chromosome names start with "chr"
    # if add_to_ins == TRUE, adds one to end of insertions, as required for VariantAnnotation

    cig_by_run <- lapply(.self$crispr_runs, function(crun) crun$cigar_labels)
    all_cigars <- unlist(cig_by_run)
    unique_cigars <- !duplicated(unlist(cig_by_run))

    co <- do.call(c, unlist(lapply(.self$crispr_runs, function(x) x$cigar_ops),
                            use.names = FALSE))[unique_cigars]
    idxs <- co != "M"

    get_gr <- function(alns){
      cigarRangesAlongReferenceSpace(cigar(alns), pos = start(alns))
    }

    ir <- do.call(c, unlist(lapply(cset$crispr_runs, function(x) get_gr(x$alns)),
                            use.names = FALSE))

    ir <- ir[unique_cigars][idxs]
    names(ir) <- all_cigars[unique_cigars]
    ir <- unlist(ir)

    if (add_to_ins){
      ins_idxs <- unlist(co[idxs] == "I")
      end(ir[ins_idxs]) <- end(ir[ins_idxs]) + 1
    }

    chrom <- as.character(seqnames(.self$target))
    if (add_chr & ! grepl('^chr', chrom)){
      chrom <- paste0("chr", chrom)
    }

    GRanges(chrom, ir)
  },

  filterVariants = function(cig_freqs = NULL, names = NULL, columns = NULL,
                            include.chimeras = TRUE){

'
Description:
  Relabels specified variants in a table of variant allele counts as 
  non-variant, e.g. variants known to exist in control samples.
  Accepts either a size, e.g. "1D", or a specific mutation, e.g. "-4:3D".
  For alleles that include one variant to be filtered and one other variant,
  the other variant will be retained.
  If SNVs are included, these will be removed entirely, but note that SNVs
  are only called in reads that do not contain an insertion/deletion variant

Input parameters:
  cig_freqs:        A table of variant allele counts 
                    (Default: NULL, i.e. .self$cigar_freqs)
  names:            Labels of variants alleles to remove (Default: NULL)
  columns:          Indices or names of control samples.  Remove all variants that 
                    occur in these columns.  (Default: NULL)
  include.chimeras: Should chimeric reads be included? (Default: TRUE)
'

    # Potential improvements:
    # Column must include only column names from .self$cigar_freqs

    warning("This function will not correctly count SNVs as variants after filtering")
    if (is.null(cig_freqs)){
      cig_freqs <- .self$.getFilteredCigarTable(include.chimeras = include.chimeras)
    }
    vars <- strsplit(rownames(cig_freqs), ",")
 
    if (length(columns) > 0){
      # Select the rownames that occur in columns
      cols <- rownames(cig_freqs)[rowSums(cig_freqs[,columns, drop=FALSE]) > 0]
    } else {
      cols <- NULL
    }

    to_remove <- c(names, cols)
    rm_snv <- grepl(.self$pars$mismatch_label, to_remove)
    has_loc <- grepl(":", to_remove) & ! rm_snv
    by_loc <- to_remove[has_loc]
    by_size <- to_remove[!(has_loc|rm_snv)]
    
    # Remove SNVs - SNVS have format SNV:-1,-5, other vars -1:3D etc
    rm_snv <- gsub(".*:(.*)", "\\1", filter.vars[rm_snv])
    temp <- gsub(".*:(.*)", "\\1", unlist(vars))
    rm_snvs <- sapply(relist(temp %in% rm_snv, vars), any)
    vars[rm_snvs] <- NULL
    cig_freqs <- cig_freqs[!rm_snvs, , drop=FALSE]

    # Reclassify indel variants as non-variant
    loc_mask <- rep(TRUE, length(unlist(vars)))
    size_mask <- loc_mask

    # Remove by location
    if (length(by_loc) > 0){
      loc_mask <- !grepl(paste(by_loc, collapse = "|"), unlist(vars))
    }
    # Remove by size - cannot be mutation combinations
    if (length(by_size) > 0){
      size_mask <- !(gsub(".*:", "", unlist(vars)) %in% by_size)
    }
    mask <- relist(loc_mask & size_mask, vars)
    vars <- as.list(IRanges::CharacterList(vars)[mask])
    vars <- lapply(vars, paste, sep = ",", collapse = ",")

    # Here: by counting as non-variant, some SNVs may be missed
    vars[vars == ""] <- .self$pars$match_label
    cig_freqs <- rowsum(cig_freqs, unlist(vars))
    cig_freqs
  },

  .getSNVs = function(min.freq = 0.25, include.chimeras = TRUE){
    # Relies on having short cigars, remove other options?

    cig_fqs <- .self$.getFilteredCigarTable(include.chimeras = include.chimeras)
    snv <- .self$pars["mismatch_label"]
    snv_nms <- rownames(cig_fqs)[grep(snv, rownames(cig_fqs))]
    all_snv_locs <- strsplit(gsub(sprintf("%s|:", snv), "", snv_nms), ",")
    snv_locs <- unique(unlist(all_snv_locs))
    asv <- unlist(all_snv_locs)
    total_count <- sum(cig_fqs)
    snv_fqs <- structure(vector(length = length(snv_locs)), names = snv_locs)
    for (snv_loc in snv_locs){
      nms <- snv_nms[sapply(relist(asv == snv_loc, all_snv_locs), any)]
      snv_count <- sum(cig_fqs[nms,])
      fq <- snv_count/total_count
      snv_fqs[snv_loc] <- fq
    }
    snv_fqs[snv_fqs >= min.freq]
  },

  mutationEfficiency = function(snv = c("non_variant", "include","exclude"),
                                include.chimeras = TRUE,
                                exclude.cols = NULL, group = NULL,
                                filter.vars = NULL, filter.cols = NULL,
                                count.alleles = FALSE, per.sample = TRUE){
'
Description:
  Calculates summary statistics for the mutation efficiency, i.e.
  the percentage of reads that contain a variant.  Reads that do not
  contain and insertion or deletion, but do contain a single nucleotide
  variant (snv) can be considered as mutated, non-mutated, or not
  included in efficiency calculations as they are ambiguous.

Input parameters:
  snv:    One of "include" (consider reads with mismatches to be mutated),
          "exclude" (do not include reads with snvs in efficiency calculations),
          and "non_variant" (consider reads with mismatches to be non-mutated).
  include.chimeras: Should chimeras be counted as variants?  (Default: TRUE)
  exclude.cols:   A list of column names to exclude from calculation, e.g. if one
                  sample is a control (default: NULL, i.e. include all columns)
  group:          A grouping variable.  Efficiency will be calculated per group,
                  instead of for individual.
                  Cannot be used with exclude.cols.
  filter.vars:    Variants that should not be counted as mutations.
  filter.cols:    Column names to be considered controls.  Variants occuring in
                  a control sample will not be counted as mutations.
  count.alleles:  If TRUE, also report statistics about the number of alleles
                  per sample/per group. (Default:  FALSE)
  per.sample:     Return efficiencies for each sample (Default: TRUE)
Return value:
  A vector of efficiency statistics per sample and overall, or a
  matrix if a group is supplied.

'
    if (! is.null(group) && ! is.null(exclude.cols)){
      stop("Only one of group or exclude.cols may be supplied")
    }

    snv <- match.arg(snv)
    freqs <- .self$cigar_freqs
    if (isTRUE(include.chimeras)){
      freqs <- .self$.getFilteredCigarTable(include.chimeras = include.chimeras)
    }

    exclude.idxs <- match(exclude.cols, colnames(freqs))
    if (any(is.na(exclude.idxs))){
      nf <- exclude.cols[is.na(exclude.idxs)]
      stop(sprintf("Column(s) %s not found in variant counts table",
                   paste(nf, collapse = " ")))
    }

    if (length(filter.vars) > 0 | length(filter.cols) > 0){
      freqs <- .self$filterVariants(cig_freqs = freqs, names = filter.vars,
                                    columns = filter.cols)
    }

    if (length(exclude.idxs) > 0){
      freqs <- freqs[,-exclude.idxs, drop = FALSE]
    }

    is_snv <- grep(.self$pars$mismatch_label, rownames(freqs))

    if (snv == "exclude"){
      if (length(is_snv) > 0) freqs <- freqs[-is_snv,,drop = FALSE]
    }

    total_seqs <- colSums(freqs)
    names(total_seqs) <- colnames(freqs)
    not_mutated <- grep(.self$pars$match_label, rownames(freqs))
    if (snv == "non_variant") not_mutated <- c(not_mutated, is_snv)

    if (length(not_mutated) > 0) freqs <- freqs[-not_mutated,,drop = FALSE]

    if (! is.null(group)){
      result <- lapply(levels(group), function(g){
        eff <- .self$.calculateEfficiency(freqs[,group == g, drop = FALSE],
                                          total_seqs, count.alleles,
                                          per.sample)
        
      })
      result <- do.call(rbind, result)
      rownames(result) <- levels(group)
      return(result)
    }

    .self$.calculateEfficiency(freqs, total_seqs, count.alleles, per.sample)
  },

  .calculateEfficiency = function(freqs, total_seqs, count.alleles, 
                                  per.sample){
     mutants <- colSums(freqs)
     ts <- total_seqs[colnames(freqs)]
     mutant_efficiency = mutants/ts * 100
     average <- mean(mutant_efficiency)
     median <- median(mutant_efficiency)
     overall <- sum(mutants)/ sum(ts) * 100
     sds <- sd(mutants)
     result <- round(c(mutant_efficiency, average, median, overall, sds, 
                       sum(ts)),2)
     names(result) <- c(colnames(freqs), "Average","Median","Overall",
                        "StDev","ReadCount")
     
     if (isTRUE(count.alleles)){
       result <- c(result, .self$.countVariantAlleles(freqs))
     }
     if (! isTRUE(per.sample)){
       result <- result[!names(result) %in% colnames(freqs)]
     }
     result
  },

  .countVariantAlleles = function(freqs){
    # how many nonzero entries per sample (absolute and percentage)
    # This step removes variants that never occur in this sample group
    freqs <- freqs[rowSums(freqs) > 0,, drop = FALSE]
    if (nrow(freqs) == 0) return(vector())
    alleles <- colSums(freqs > 0)
    alleles_pc <- alleles/nrow(freqs)*100
    
    # Rate of new alleles (how many sequences per allele)
    # Produces "NA" if any samples have 0 alleles
    allele_rate <- colSums(freqs)/alleles
    nvars <- nrow(freqs)
    
    result <- c(mean(alleles), mean(alleles_pc), mean(allele_rate), min(alleles),
                max(alleles), nvars)
    result <- sapply(result, round, 2)
    names(result) <- c("AvAlleles","AvPctAlleles","SeqsPerAllele", 
                       "MinAlleles", "MaxAlleles","TotalAlleles")
    result
  },

  .classifyVariantsByType = function(){
'
Description:
  Classifies variants as insertions, deletions, or complex (combinations).
  In development
'
    
    # Classifies variants as reference, mismatch, insertion, deletion
    # or insertion+deletion
    vars <- rep(NA, nrow(.self$cigar_freqs))
    is_snv <- grepl(.self$pars$mismatch_label, rownames(.self$cigar_freqs))
    is_ref <- grepl(.self$pars$match_label, rownames(.self$cigar_freqs))
    is_ins <- grepl("I", rownames(.self$cigar_freqs))
    is_del <- grepl("D", rownames(.self$cigar_freqs))
    ins_and_del <- is_ins & is_del
    vars[is_ref] <- .self$pars$match_label
    vars[is_snv] <- .self$pars$mismatch_label
    vars[is_ins] <- "insertion"
    vars[is_del] <- "deletion"
    vars[ins_and_del] <- "insertion/deletion"
    vars
  },

  classifyVariantsByLoc = function(txdb, add_chr = TRUE, verbose = TRUE){
  '
Description:
  Uses the VariantAnnotation package to look up the location of the
  variants.  VariantAnnotation allows multiple classification tags per variant,
  this function returns a single tag.  The following preference order is used:
  spliceSite > coding > intron > fiveUTR > threeUTR > promoter > intergenic

Input parameters:
  txdb:     A BSgenome transcription database
  add_chr:  Add "chr" to chromosome names to make compatible with UCSC (default: TRUE)
  verbose:  Print progress (default: TRUE)

Return value:
  A vector of classification tags, matching the rownames of .self$cigar_freqs
  (the variant count table)
  '

    if (verbose) message("Looking up variant locations\n")

    stopifnot(requireNamespace("VariantAnnotation"))

    gr <- .self$.getUniqueIndelRanges(add_chr)
    locs <- VariantAnnotation::locateVariants(gr, txdb, VariantAnnotation::AllVariants())
    if (isTRUE(verbose)) message("Classifying variants\n")

    locs_codes <- paste(seqnames(locs), start(locs), end(locs), sep = "_")
    # Note that all indels have the same range
    indel_codes <- paste(seqnames(gr), start(gr), end(gr), sep = "_")
    indel_to_loc <- lapply(indel_codes, function(x) locs$LOCATION[which(locs_codes == x)])

    var_levels <- c("spliceSite","coding","intron","fiveUTR",
                    "threeUTR","promoter", "intergenic")
    result <- unlist(lapply(indel_to_loc, function(x){
      y <- factor(x,levels = var_levels)
      var_levels[min(as.numeric(y))]}))
    names(result) <- names(gr)

    classification <- rep("", nrow(.self$cigar_freqs))
    no_var <- grep(.self$pars$match_label, rownames(.self$cigar_freqs))
    classification[no_var] <- .self$pars$match_label
    snv <- grep(.self$pars$mismatch_label, rownames(.self$cigar_freqs))
    classification[snv] <- .self$pars$mismatch_label

    ord <- match(names(result), rownames(.self$cigar_freqs))
    classification[ord] <- result
    names(classification) <- rownames(.self$cigar_freqs)

    classification
  },

  classifyCodingBySize = function(var_type, cutoff = 10){
'
Description:
    This is a naive classification of variants as frameshift or in-frame
    Coding indels are summed, and indels with sum divisible by 3 are
    considered frameshift.  Note that this may not be correct for variants
    that span an intron-exon boundary
Input paramters:
    var_type:   A vector of var_type. Only variants with var_type == "coding" 
                are considered.  Intended to work with classifyVariantsByLoc
    cutoff:     Variants are divided into those less than and greater 
                than "cutoff" (Default: 10)
Result:
    A character vector with a classification for each variant allele
'

    is_coding <- var_type == "coding"

    indels <- .self$cigar_freqs[is_coding,,drop = FALSE]
    if (length(indels) > 0){

      temp <- lapply(rownames(indels), function(x) strsplit(x, ",")[[1]])
      indel_grp <- rep(c(1:nrow(indels)), lapply(temp, length))
      indel_ln <- rowsum(as.numeric(gsub("^.*:([0-9]+)[DI]", "\\1", unlist(temp))), indel_grp)

      inframe <- indel_ln %% 3 == 0
      is_short <- indel_ln < cutoff

      indel_grp <- rep(sprintf("inframe indel < %s", cutoff), nrow(indels))
      indel_grp[is_short &! inframe] <- sprintf("frameshift indel < %s", cutoff)
      indel_grp[!is_short & inframe] <- sprintf("inframe indel > %s", cutoff)
      indel_grp[!is_short & !inframe] <- sprintf("frameshift indel > %s", cutoff)
      var_type[is_coding] <- indel_grp
    }

    var_type
  },

  heatmapCigarFreqs = function(as.percent = TRUE, x.size = 8, y.size = 8,
                               x.axis.title = NULL, x.angle = 90,
                               min.freq = 0, min.count = 0,
                               top.n = nrow(.self$cigar_freqs),
                               type = c("counts", "proportions"), ...){

    '
Description:
    Internal method for CrispRVariants:plotFreqHeatmap, optionally filters the table
    of variants, then a table of variant counts, coloured by counts or proportions.

Input parameters:
    as.percent:   Should colours represent the percentage of reads per sample 
                  (TRUE) or the actual counts (FALSE)?  (Default: TRUE)
    x.size:       Font size for x axis labels (Default: 8)
    y.size:       Font size for y axis labels (Default: 8)
    x.axis.title: Title for x axis
    min.freq:     Include only variants with frequency at least min.freq in at 
                  least one sample
    min.count:    Include only variants with count at least min.count in at 
                  least one sample
    top.n:        Include only the n most common variants
    type:         Should labels show counts or proportions?  (Default: counts)
    ...:

Return value:
    A ggplot2 plot object.  Call "print(obj)" to display

See also:
    CrispRVariants::plotFreqHeatmap
    '

    # Doesn't currently allow option to exclude non-variant

    cig_freqs <- .self$.getFilteredCigarTable(top.n, min.count, min.freq,
                                              type = type)
    header <- NULL
    type <- match.arg(type)
    if (type == "counts"){ header <- colSums(.self$.getFilteredCigarTable())
    } else if (type == "proportions"){
      header <- round(colSums(cig_freqs), 2)
      cig_freqs <- round(cig_freqs, 2)
    }

    p <- plotFreqHeatmap(cig_freqs, header = header, as.percent = as.percent,
                         x.size = x.size, y.size = y.size,
                         x.axis.title = x.axis.title,
                         x.angle = x.angle, add.other = TRUE, ...)
    return(p)
  },

  plotVariants = function(min.freq = 0, min.count = 0, top.n = nrow(.self$cigar_freqs),
                   renumbered = .self$pars["renumbered"], add.other = add.other, ...){
'
Description:
  Internal method for CrispRVariants:plotAlignments, optionally filters the table
  of variants, then plots variants with respect to the reference sequence,
  collapsing insertions and displaying insertion sequences below the plot.

Input parameters:
  min.freq:         i(%) include variants that occur in at least i% of reads
                    in at least one sample
  min.count         i (integer) include variants that occur at leas i times in
                    at least one sample
  top.n:            n (integer) Plot only the n most frequent variants
                    (default: plot all)
                    Note that if there are ties in variant ranks,
                    top.n only includes ties with all members ranking <= top.n
  renumbered:       If TRUE, the x-axis is numbered with respect to the target
                    (cut) site.  If FALSE, x-axis shows genomic locations.
                    (default: TRUE)
  add.other         Add a blank row named "Other" for chimeric alignments,
                    if there are any (Default: TRUE)
  ...               additional arguments for plotAlignments

Return value:
  A ggplot2 plot object.  Call "print(obj)" to display
'

    cig_freqs <- .self$.getFilteredCigarTable(top.n, min.count, min.freq)
    # If there are no chimeric alignments, drop "Other"
    if ("Other" %in% rownames(cig_freqs)){
      cig_freqs <- cig_freqs[rownames(cig_freqs) != "Other",, drop = FALSE]
    } else {
      add.other <- FALSE
    }

    alns <- .self$makePairwiseAlns(cig_freqs)
    dots <- list(...)

    tloc <- ifelse(is.na(.self$pars$target.loc), 17, .self$pars$target.loc)
    ins.sites <- .self$insertion_sites

    # Consistency - is it still possible to print w.r.t the forward strand?
    # If the strand is -ve, the region will be reverse complemented,
    # "start" for the plot must be reversed
    if (as.character(strand(.self$target) == "-")){
      temp <- seq_along(1:width(.self$target))
      starts <- rev(temp)
      names(starts) <- temp
      ins.sites$start <- starts[ins.sites$start]
    }

    if (isTRUE(renumbered)){
      genomic_coords <- c(start(.self$target):end(.self$target))
      target_coords <- .self$genome_to_target[as.character(genomic_coords)]
      if (as.character(strand(.self$target)) == "-"){
        target_coords <- rev(target_coords)
      }
      xbreaks = which(target_coords %% 5 == 0 | abs(target_coords) == 1)
      target_coords <- target_coords[xbreaks]

      args <- list(obj = .self$ref, alns = alns, ins.sites = ins.sites,
                   xtick.labs = target_coords, xtick.breaks = xbreaks,
                   target.loc =  tloc, add.other = add.other)
    } else {
      args <- list(obj = .self$ref, alns = alns, ins.sites = ins.sites,
                   target.loc = tloc, add.other = add.other)
    }

    args <- modifyList(args, dots)
    p <- do.call(plotAlignments, args)

    return(p)
  },

  .getInsertions = function(with_cigars = TRUE){
    # Used by plotVariants for getting a table of insertions

    if (with_cigars == FALSE){
      all_ins <- do.call(rbind, lapply(.self$crispr_runs, function(x) x$insertions))
    } else {
      all_ins <- do.call(rbind, lapply(.self$crispr_runs, function(x) {
        ik <- x$ins_key
        v <- data.frame(ik, x$cigar_labels[as.integer(names(ik))])
        v <- v[!duplicated(v),]
        v <- v[order(v$ik),]
        cbind(x$insertions[v[,1],], cigar = v[,2])
      }))
    }
    if (nrow(all_ins) == 0) {
      .self$field("insertion_sites", all_ins)
      return()
    }
    new_ins_sites <- all_ins[order(all_ins$start, all_ins$seq),, drop = FALSE]
    .self$field("insertion_sites", new_ins_sites)
  },

  makePairwiseAlns = function(cig_freqs = .self$cigar_freqs, ...){
'
Description:
Get variants by their cigar string, make the pairwise alignments for the consensus
sequence for each variant allele

Input parameters:
cig_freqs:  A table of variant allele frequencies (by default: .self$cigar_freqs,
            but could also be filtered)
...:        Extra arguments for CrispRVariants::seqsToAln, e.g. which symbol
            should be used for representing deleted bases
'

    # The short cigars (not renumbered) do not have enough information,
    # use the full cigars for sorting
    # Do this just for the alns to be displayed?

    cigs <- unlist(lapply(.self$crispr_runs, function(x) cigar(x$alns)), use.names = FALSE)
    cig_labels <- unlist(lapply(.self$crispr_runs, function(x) x$cigar_labels), use.names = FALSE)

    names(cigs) <- cig_labels # calling by name with duplicates returns the first match

    splits <- split(seq_along(cig_labels), cig_labels)
    splits <- splits[match(rownames(cig_freqs), names(splits))]

    splits_labels <- names(splits)
    names(splits) <- cigs[names(splits)]
    all_d <- grep("M", names(splits), invert = TRUE)

    x <- lapply(.self$crispr_runs, function(x) x$alns)
    all_alns <- do.call(c, unlist(x, use.names = FALSE))

    seqs <- c()
    starts <- c()

    for (i in seq_along(splits)){
      if (i %in% all_d){
        seqs[i] <- ""
      } else {
        idxs <- splits[[i]]
        cig <- names(splits[[i]])
        seqs[i] <- consensusString(mcols(all_alns[idxs])$seq)
      }
      start <- unique(start(all_alns[idxs]))
      if (length(start) > 1)
        stop("Sequences with the same cigar string have different starting locations.
             This case is not implemented yet.")
      starts[i] <- start[1]
    }

    seqs <- Biostrings::DNAStringSet(seqs)
    alns <- seqsToAln(names(splits), seqs, target = .self$target,
                 aln_start = starts, ...)

    names(alns) <- splits_labels
    alns
  },

  .genomeToTargetLocs = function(target.loc, target_start, target_end, rc = FALSE){
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

    gs <- min(sapply(.self$crispr_runs, function(x) min(start(x$alns))))
    ge <- max(sapply(.self$crispr_runs, function(x) max(end(x$alns))))

    if (isTRUE(rc)){
      tg <- target_end - (target.loc - 1)
      new_numbering <- rev(c(seq(-1*(ge - (tg -1)),-1), c(1:(tg - gs))))
      names(new_numbering) <- c(gs:ge)

    } else {
      tg <- target_start + target.loc - 1
      new_numbering <- c(seq(-1*(tg - (gs-1)),-1), c(1:(ge - tg)))
      names(new_numbering) <- c(gs:ge)
    }
    .self$field("genome_to_target", new_numbering)
    new_numbering
  }
)