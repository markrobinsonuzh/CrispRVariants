#'@title Find chimeric reads
#'@description Find chimeric reads, assuming that the GAlignments object
#' does not contain multimapping reads. That is, read names that appear
#' more than ones in the file are considered chimeras.  Chimeric reads
#' are reads that cannot be mapped as a single, linear alignment.  Reads 
#' from structual rearrangements such as inversions can be mapped as chimeras.
#' Note that the indices of all chimeric reads are returned, these are not 
#' separated into individual chimeric sets.
#'@param bam  A GAlignments object, must include names
#'@param by.flag Can the chimeras be detected just using the supplementary
#'alignment flag?  (Default: FALSE).  If TRUE, detects supplementary alignments 
#'and returns reads with the same name as a supplementary alignment (quicker).
#'If FALSE, all alignments with duplicated names are returned.
#'@author Helen Lindsay
#'@return A vector of indices of chimeric sequences within the original bam
#'@seealso \code{\link{plotChimeras}} for plotting chimeric alignment sets.
#'@export
#'@examples
#'bam_fname <- system.file("extdata", "gol_F1_clutch_2_embryo_4_s.bam",
#'                          package = "crispRvariants")
#'bam <- GenomicAlignments::readGAlignments(bam_fname, use.names = TRUE)
#'chimera_indices <- findChimeras(bam)
#'chimeras <- bam[chimera_indices]
findChimeras <- function(bam, by.flag = FALSE){
  if (by.flag == TRUE & ! "flag" %in% names(mcols(bam))){
    stop("If 'by.cigar' is TRUE, bam must have a metadata column 'flag'")
  }
  if (by.flag == TRUE){
    suppl <- bitwAnd(mcols(bam)$flag, 2048)
    ch_names <- names(bam[suppl])
    chimera_idxs <- names(bam) %in% ch_names
  } else {
    chimera_idxs <- which((duplicated(names(bam)) | 
                           duplicated(names(bam), fromLast = TRUE)))
  }
  chimera_idxs <- chimera_idxs[order(as.factor(names(bam)[chimera_idxs]))] 
  return(chimera_idxs)
}

mergeChimeras <- function(bam){
  # Note: this does not check whether chimeras are mergeable (strands, chromosomes) - do not use!
  # Cannot incorporate actual insertions, just the unmapped length
  
  # To do: store the chimeric segments that aren't on the target chromosome /
  # don't end up in the target region
  
  chimera_idxs <- findChimeras(bam)
  chimeras <- bam[chimera_idxs]
  nms <- rle(names(bam)[chimera_idxs]) 
  nch <- length(chimera_idxs)
  
  # Find start and end points of chimeras within the set
  ch_ends <- cumsum(nms$lengths)
  change_pts <- c(1, ch_ends + 1)[1:length(ch_ends)]
  
  cigars <- cigar(bam)[chimera_idxs]
  unclipped <- gsub("[0-9]+[HS]$", "", gsub("^[0-9]+[HS]", "", cigars))
  first_aligned <- rep(1, length(cigars))
  clipped_start <- grepl("^[0-9]+[HS]", cigars)    
  # Note: first aligned refers to the original read, not the clipped read
  first_aligned[clipped_start] <- as.numeric(gsub("[HS].*", "", cigars[clipped_start])) + 1 
  cig_ranges <- cigarRangesAlongQuerySpace(unclipped)
  last_aligned <- sum(width(cig_ranges))  
  last_aligned <- last_aligned  + first_aligned - 1
  
  # How much do reads overlap?  Cut bases that align to two segs from one seg
  read_gaps <- c(0,first_aligned[-1] -  last_aligned[-length(last_aligned)] - 1)
  read_gaps[change_pts] <- 0
  to_cut <- read_gaps < 0
  
  new_cigars <- unclipped
  # For first member of chimera: get everything up to and including the last M
  new_cigars[change_pts] <- gsub("(^.*M)[0-9]+[HS]","\\1", cigars[change_pts])
  # For last member of chimera: get everything except for clipping at the start
  new_cigars[ch_ends] <- gsub("H","S",gsub("^[0-9]+[HS](.*)", "\\1", cigars[ch_ends]))
  
  # Get the first alignment operation for alignments that will be trimmed
  first_range <- as.numeric(gsub('M.*', "", new_cigars[to_cut]))
  # To do: (?)
  # check that value to cut isn't longer that the first operation
  
  # Cut overlap off rightmost read, adjust genomic coordinates
  first_range <- first_range + read_gaps[to_cut]
  new_cigars[to_cut] <- paste0(first_range, gsub('[0-9]+(M.*)', "\\1", new_cigars[to_cut]))
  
  # Find new genomic start coordinates, accounting for clipping in original read
  genomic_starts <- start(chimeras)
  genomic_starts[to_cut] <- genomic_starts[to_cut] - read_gaps[to_cut]
  
  # Stick cigars together padding genomic gaps with deletions
  ggaps <- c(sprintf("%sD", (start(chimeras[-1]) - end(chimeras[-length(chimeras)])) +1),0)
  ggaps[ch_ends] <- ""
  ggaps[ggaps == "0D"] <- ""
  
  # add insertions for read gaps
  rgaps <- sprintf("%sI", read_gaps)
  rgaps[change_pts] <- ""
  rgaps[rgaps == "0I"] <- ""
  
  new_cigars <- paste0(rgaps, new_cigars)
  new_cigars <- paste0(new_cigars, ggaps)

  # Paste all cigars from a chimeric set together
  new_cigars <- aggregate(new_cigars, list(names(chimeras)), FUN = paste0, collapse = "")$x
  seqs <- mcols(chimeras[change_pts])$seq
  ends <- aggregate(end(chimeras), list(names(chimeras)), FUN = max)$x
  chim <- GRanges(seqnames(chimeras)[change_pts], 
                  IRanges(start(chimeras)[change_pts], ends), 
                  cigar = new_cigars, strand = strand(chimeras[change_pts]))
  mcols(chim)$seq <- seqs
  mcols(chim)$flag <- mcols(chimeras)$flag[change_pts]
  names(chim) <- names(chimeras)[change_pts]
  chim <- as(chim, "GAlignments")
  nchim <- bam[setdiff(seq_along(bam), chimera_idxs)]
  result <- c(chim, nchim)
  return(result)
  
}