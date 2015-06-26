#'@title Display a dot plot of chimeric alignments
#'
#'@description Produces a dot plot of a set of chimeric alignments.  For chimeric alignments,
#'a single read is split into several, possibly overlapping alignmed blocks.
#'Aligned sections of chimeric reads can be separated by large genomic distances, 
#'or on separate chromosomes.  plotChimeras produces a dot plot, each aligned block 
#'highlighted, and chromosomes shown in different colours. Large gaps between 
#'aligned segments are collapsed and indicated on the plot with horizontal lines.  
#'The X-axis shows each base of the entire read. Note that the mapping to the fwd strand 
#'is shown if all strands agree.  The chimeric alignments must be sorted!   
#'
#'@param chimeric.alns A GAlignments object containing only the chimeric
#' reads to be plotted 
#'@param max.gap If aligned segments are separated by more than max.gap, 
# with respect to the genome, a gap in the y-axis will be introduced (Default 10)
#'@param tick.sep How many bases should separate tick labels on plot.  Default 20.
#'@param text.size Size of X and Y tick labels on plot.  Default 12
#'@param title.size Size of X and Y axis labels on plot.  Default 16
#'@param gap.pad How much should aligned blocks be separated by?  (Default: 20)
#'@param legend.title Title for the legend.  Default "Chromosome"
#'@param xangle Angle for x axis text (Default 90, i.e vertical)
#'@param wrt.forward Should chimeric alignments where all members map to the 
#'negative strand be displayed with respect to the forward strand, i.e. as the 
#'cigar strand is written (TRUE), or the negative strand (FALSE) (Default: FALSE)
#'@param annotations A list of GRanges.  Any that overlap with the chimeric alignments
#'are highlighed in the plot.  
#'@param annotate.within annot_aln ranges in "annotations" within n bases of a chimeric 
#'alignment (Default 50)
#'@return A ggplot2 dotplot of the chimeric alignments versus the reference sequence 
#'@seealso \code{\link{findChimeras}} for finding chimeric alignment sets.
#'@author Helen Lindsay
#'@export
#'@examples
#'bam_fname <- system.file("extdata", "gol_F1_clutch_2_embryo_4_s.bam",
#'                          package = "crispRvariants")
#'bam <- GenomicAlignments::readGAlignments(bam_fname, use.names = TRUE)
#'# Choose a single chimeric read set to plot:
#'chimeras <- bam[names(bam) == "AB3092"]
#'
#'# This read aligns in 3 pieces, all on chromosome 18.  
#'# The plot shows the alignment annot_alns a small duplication and 
#'# a long gap.
#'plotChimeras(chimeras)
plotChimeras <- function(chimeric.alns, max.gap = 10, tick.sep = 20, 
                         text.size = 10,  title.size = 16, gap.pad = 20,
                         legend.title = "Chromosome", xangle = 90, 
                         wrt.forward = FALSE, annotate.within = 20,
                         annotations = GenomicRanges::GRanges()){
  # TO DO?
  # - extend y-axis to annotate regions that are close but not spanned
  # - label the chromosomal regions
  # Count soft-clipped bases at ends?
  
  # Sort chromosomes by minimum read start of any segment on each chr
  temp_cigs <- cigarRangesAlongQuerySpace(cigar(chimeric.alns))
  is_match <- CharacterList(explodeCigarOps(cigar(chimeric.alns))) == "M"
  sq_ord <- seqnames(chimeric.alns[order( min(start(temp_cigs[is_match])))])
  seqlevels(chimeric.alns) <- unique(as.character(sq_ord@values))
  chimeric.alns <- chimeric.alns[order(seqnames(chimeric.alns))]
    
  # Get M ranges wrt genome and read
  cigars <- cigar(chimeric.alns)
  genomic_locs <- as(chimeric.alns, "GRanges")
  is_plus <- as.vector(strand(genomic_locs) == "+")
  two_strands <- length(unique(is_plus)) > 1
 
  # For reference ranges, shift to actual genomic starting locations
  ref_ranges <- cigarRangesAlongReferenceSpace(cigars)
  ref_ranges <- shift(ref_ranges, start(genomic_locs) -1)
  
  # Entirely negative strand chimeras may be displayed as aligned to reference
  # (wrt.forward = TRUE).  Default is wrt start of read
  if (two_strands == FALSE & wrt.forward == FALSE & ! any(is_plus)){
    cigars <- reverseCigarv(cigars)
    ref_ranges <- relist(rev(unlist(rev(ref_ranges))), ref_ranges)
  }
  
  # Extract match ranges, these will be plotted
  ops <- CharacterList(explodeCigarOps(cigars))
  query_ranges <- cigarRangesAlongQuerySpace(cigars)
  
  # Find all "M" operations (runs of aligned bases)
  mm <- ops == "M"
  mm_idxs <- rep(1:length(genomic_locs), sum(mm))
  m_ref <- ref_ranges[mm]
  m_qry <- query_ranges[mm]
  
  # Split by strand and add offset for hard clipping
  hclipl <- as.numeric(gsub(".*[A-Z].*", 0, gsub("[H].*","", cigars)))
  m_qry[is_plus] <- GenomicRanges::shift(m_qry[is_plus], hclipl[is_plus])
  
  # If the alignment annot_alns segments to both strands, display wrt ref (+)
  if (two_strands) { 
    # Want to get all query ranges with respect to the query on the forward strand
    # (As is, cigar ranges are with respect to the genome)
    hsclipr <- as.numeric(gsub("^$", 0, gsub('.*[M]|[HS]$', "", cigars))) 
    
    # Find the distance from the end of the -ve strand M ranges to the start of the last op
    dist_to_last <- (max(end(query_ranges)[mm]) - end(query_ranges))[mm][!is_plus]
    
    # Remove the offset from the left, add the offset from the right
    m_qry[!is_plus] <- GenomicRanges::shift(m_qry[!is_plus], (-1*start(m_qry[!is_plus])+1))
    m_qry[!is_plus] <- GenomicRanges::shift(m_qry[!is_plus], hsclipr[!is_plus] + dist_to_last) 
    
  } else {
  # Else if only -ve, display wrt negative
    m_qry[!is_plus] <- GenomicRanges::shift(m_qry[!is_plus], hclipl[!is_plus])
  }
  
  # Ord is the order of the "M" segments (can be more than number of alignments)
  plus <- rep(which(is_plus), lapply(m_qry[is_plus], length))
  minus <- rep(which(!is_plus), lapply(m_qry[!is_plus], length))
  ord <- order(c(plus,minus))

  xs <- mapply(seq, unlist(start(m_qry[is_plus])), unlist(end(m_qry[is_plus])), 
               SIMPLIFY = FALSE)

  if (wrt.forward == FALSE){
    xs <- c(xs, mapply(function(x,y) rev(seq(x,y)), unlist(start(m_qry[!is_plus])), 
          unlist(end(m_qry[!is_plus])), SIMPLIFY = FALSE))[ord]
  } else {
    xs <- c(xs, mapply(function(x,y) seq(x,y), unlist(start(m_qry[!is_plus])), 
          unlist(end(m_qry[!is_plus])), SIMPLIFY = FALSE))
  }
  
  ys <- mapply(seq, unlist(start(m_ref)), unlist(end(m_ref)), SIMPLIFY = FALSE)
  chrs <- seqnames(genomic_locs)
  
  #____________________________
  # Introduce gaps for aligned segments separated by more than max.gap
  # NOTE: assume chimeric alns are sorted (ys segment cannot be before previous)
  
  y_lns <- sapply(ys, length)
  chr_to_ranges <- as.character(seqnames(genomic_locs)[mm_idxs])
  m_granges <- GRanges(chr_to_ranges, unlist(m_ref))
  y_blocks <- reduce(m_granges, min.gapwidth = max.gap)
  # Blocks that are not merged should have a gap added
  
  # SIMPLIFY OFFSETS BY CALCULATING WRT BLOCK
  
  offsets <- rep(0, length(m_granges))
  ys_to_block <- subjectHits(findOverlaps(m_granges, y_blocks))
  offsets[!duplicated(ys_to_block) & ys_to_block!= 1] <- gap.pad
  # map offsets to y_block
  offsets <- offsets[!duplicated(ys_to_block)]
  ycoords <- 1:sum(width(y_blocks))
  names(ycoords) <- unlist(mapply(seq, start(y_blocks), end(y_blocks)))
  # map y_block offsets to correct y segment, add to y coordinates
  offsets <- rep(cumsum(offsets)[ys_to_block], y_lns)
  ycoords <- ycoords[as.character(unlist(ys))] + offsets
  
  #________________________________
  # Make boxes around each segment
  
  y_sums <- cumsum(y_lns)
  n_segs <- length(unlist(m_qry))
  
  # Setup data and plot
  pt_coords <- data.frame(x = unlist(xs), ylabs = unlist(ys))
  pt_coords$y <- ycoords 
  
  ymax <- pt_coords$y[y_sums]
  ymin <- na.omit(c(1, pt_coords$y[y_sums +1]))
  
  # Make boxes around single cigar segments to distinguish indels from chimeras
  m_qry_range <- range(m_qry)
  m_qry_start <- unlist(start(m_qry_range))
  m_qry_end <- unlist(end(m_qry_range))
  aln_to_ranges <- relist(1:n_segs, m_qry)
  
  # This requires the ys to be in the original order
  pt_ys <- relist(pt_coords$y, ys)
  box_ranges <-  t(sapply(1:length(m_qry), function(i){
      aln_ys <- unlist(pt_ys[aln_to_ranges[[i]]])
      c(min(aln_ys), max(aln_ys))
  }))
  
  box_ymins <- box_ranges[,1]
  box_ymaxs <- box_ranges[,2]
  
  box_coords <- data.frame("xmin" = m_qry_start, "xmax" = m_qry_end, 
                           "chrs" = seqnames(genomic_locs),
                           "ymin"  = box_ymins, "ymax" = box_ymaxs)
  
  #____________________________
  # Get coordinates of large breaks for shading
  
  nblocks <- length(y_blocks)
  if (nblocks > 1){
    gap_starts <- end(y_blocks)[1:(nblocks -1)]
    gap_ends <- start(y_blocks)[2:nblocks]
    chr_box_coords <- data.frame(ymin = ycoords[as.character(gap_starts)], 
                                 ymax = ycoords[as.character(gap_ends)])
  } else {
    chr_box_coords <- data.frame()
  }
  
  #____________________________
  # Make tick labels and locations
  
  int_ys <- as.integer(names(ycoords))
  tick_labs <- int_ys[int_ys %% tick.sep == 0]

  # Labels for breakpoints are added if they are not too close to existing 
  # axis labels.  Not elegant, but prevents overlapping labels
  chr_labs <- c(start(y_blocks), end(y_blocks))
  min_dists <- sapply(chr_labs, function(x) min(abs(x-tick_labs)))
  chr_labs <- chr_labs[min_dists > 2]
                 
  tick_labs <- unique(c(chr_labs, tick_labs))
  tick_labs <- tick_labs[order(tick_labs)]
  tick_pos <- ycoords[as.character(tick_labs)]

  xbreaks <- seq(min(pt_coords$x), max(pt_coords$x), by = tick.sep)
 
  #____________________________
  # Plotting
  
  p <- ggplot(pt_coords, aes(x=x,y=y)) + geom_point() +
       geom_rect(data = box_coords, aes(xmin = xmin, xmax = xmax, 
                  ymin = ymin, ymax = ymax, x = NULL, y = NULL, fill = chrs,
                  colour = chrs), alpha = 0.25) + 
       xlab("Read location") + ylab("Chromosomal location") +
       guides(fill = guide_legend(title = legend.title), 
              colour = guide_legend(title = legend.title)) +
       theme_bw() + theme(axis.title.y=element_text(vjust = 2, size = title.size), 
                          axis.title.x=element_text(vjust = -1, size = title.size),
                          axis.text.y=element_text(size=text.size),
                          axis.text.x=element_text(size=text.size, angle = xangle),
                          plot.margin = grid::unit(c(1, 1, 1, 1), "lines")) 
  if (nrow(chr_box_coords) > 0){                         
    p <- p + geom_rect(data = chr_box_coords, aes(ymin = ymin, 
                       ymax = ymax, y = NULL, x = NULL), 
                       xmin = min(pt_coords$x), xmax = max(pt_coords$x), 
                       fill = "gray", alpha = 0.2, colour = "gray", linetype = "dotted") +
         scale_x_continuous(expand = c(0,0), breaks = xbreaks) 
  }
  
  #____________________________
  # If a list of points to annotate is provided, annote overlapping points  
  
  # Annotate the aligned ranges
  aln_ranges <- GRanges(chr_to_ranges, unlist(m_ref))
  annot_aln <- annotations[queryHits(findOverlaps(annotations, aln_ranges,
                                  type = "within"))]
  
  gap_ranges <- gaps(aln_ranges, start = min(start(aln_ranges)))
  a_to_gap <- findOverlaps(annotations, gap_ranges)
  
  # How far is annotation through gap range as a percentage 
  a_starts <- start(annotations[queryHits(a_to_gap)])
  g_starts <- start(gap_ranges[subjectHits(a_to_gap)])
  g_to_y <- ycoords[as.character(g_starts - 1)]
  seps <- width(gap_ranges[subjectHits(a_to_gap)])
  seps[seps >= max.gap] <- tick.sep 
  g_to_y <- g_to_y + ((a_starts - g_starts) / 
                        width(gap_ranges[subjectHits(a_to_gap)])) * seps + 1
  g_to_y_nms <- mcols(annotations[queryHits(a_to_gap)])$name
  
  #Annotate boundaries between chromosomal blocks
  lblocks <- flank(y_blocks, annotate.within)
  lblocks <- setdiff(disjoin(c(lblocks, y_blocks, gap_ranges)), c(y_blocks, gap_ranges))
  annot_left <- findOverlaps(annotations, lblocks, type = "within")
  rblocks <- flank(y_blocks, annotate.within, start = FALSE)
  rblocks <- setdiff(disjoin(c(rblocks, y_blocks, gap_ranges)), c(y_blocks, gap_ranges))
  annot_right <- findOverlaps(annotations, rblocks, type = "within")
  
  # Distance isn't very meaningful in gaps between chromosomes,
  # so just place the annotation on the correct side
  nearest_left <- end(lblocks[subjectHits(annot_left)]) + 1
  left_locs <- ycoords[as.character(nearest_left)] - 
          min((nearest_left - start(annotations[queryHits(annot_left)])), tick.sep/4) 
  
  nearest_right <- start(rblocks[subjectHits(annot_right)]) - 1
  right_locs <- ycoords[as.character(nearest_right)] + 
        min((start(annotations[queryHits(annot_right)]) - nearest_right), tick.sep/4) 

  flank_nms <- mcols(annotations[c(queryHits(annot_left), queryHits(annot_right))])$name
  others <- c(g_to_y, left_locs, right_locs)
  
  if (length(annot_aln) > 0 | length(others) > 0){ 
    
    annot <- data.frame(yint = c(ycoords[as.character(start(annot_aln))],
                                 others))
    
    p <- p + geom_hline(data = annot, aes(yintercept = yint), linetype = "longdash", 
                        color = "red", size = 0.75)
    if ("name" %in% names(mcols(annotations))){
      
      vjust <- ifelse(annot$yint >= max(ycoords), 1,0)
      annot_lab <- data.frame(label = c(annot_aln$name, g_to_y_nms, flank_nms), 
                              vjust = vjust,
                              yint = annot$yint)
      p <- p + geom_text(data = annot_lab, x = min(pt_coords$x) + 1, 
                         aes(label = label, y = yint, vjust = vjust),
                         color = "red", hjust = 0)
    }
  }
  #____________________________
  
  ylower <- min(0, c(ycoords[as.character(start(annot_aln))], left_locs-1))
  yupper <- max(max(ymax) + 1, right_locs)
  
  p <- p + scale_y_continuous(expand = c(0,0), breaks = tick_pos, labels = tick_labs,
                              limits = c(ylower, yupper))
  p
}