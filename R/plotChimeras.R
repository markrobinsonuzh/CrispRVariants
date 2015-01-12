#'Display a dot plot of chimeric alignments
#'
#'Produces a dot plot of a set of chimeric alignments.  For chimeric alignments,
#'a single read is split into several, possibly overlapping alignmed blocks.
#'Aligned sections of chimeric reads can be separated by large genomic distances, 
#'or on separate chromosomes.  plotChimeras produces a dot plot, each aligned block 
#'highlighted, and chromosomes shown in different colours. Large gaps between 
#'aligned segments are collapsed and indicated on the plot with horizontal lines.  
#'The X-axis shows each base of the entire read.    
#'
#'@param chimeric_alns A GAlignments object containing only the chimeric
#' reads to be plotted 
#'@param max_gap The maximum gap between aligned blocks allowed 
#'aligned before blocks are collapsed, in bases.  Default 10
#'@param tick_sep How many bases should separate tick labels on plot.  Default 20.
#'@param text_size Size of X and Y tick labels on plot.  Default 12
#'@param title_size Size of X and Y axis labels on plot.  Default 16
#'@param legend_title Title for the legend.  Default "Chromosome"
#'@author Helen Lindsay
#'
#'@import GenomicAlignments
#'@import ggplot2
#'@import grid
#'@import gridExtra
#'@import GenomicRanges
#'@import IRanges
#'@import methods
#'@import parallel
#'@importFrom reshape2 melt
#'@export
#'
plotChimeras <- function(chimeric_alns, max_gap = 10, tick_sep = 20, text_size = 12,  
                         title_size = 16, legend_title = "Chromosome"){
  # max_gap: if genomic segments are separated by more than max_gap, 
  #          a gap in the y-axis will be introduced
  
  # To do: allow a region of interest to be annotated
  # note that the mapping to the fwd strand is shown if all strands agree
  cigars <- cigar(chimeric_alns)
  genomic_locs <- as(chimeric_alns, "GRanges")
  
  ops <- CharacterList(explodeCigarOps(cigars))
  query_ranges <- cigarRangesAlongQuerySpace(cigars)
  ref_ranges <- cigarRangesAlongReferenceSpace(cigars)
  ref_ranges <- shift(ref_ranges, start(genomic_locs) -1)
  mm <- ops == "M"
  mm_idxs <- rep(1:length(genomic_locs), sum(mm))
  m_ref <- ref_ranges[mm]
  m_qry <- query_ranges[mm]
  
  # Split by strand and add offset for hard clipping
  hclipl <- as.numeric(gsub(".*[A-Z].*", 0, gsub("[H].*","", cigars)))
  is_plus <- as.vector(strand(genomic_locs) == "+")
  two_strands <- length(unique(is_plus)) > 1

  m_qry[is_plus] <- shift(m_qry[is_plus], hclipl[is_plus])
  
  if (two_strands) { 
    hclipr <- as.numeric(gsub("^$", 0, gsub('.*[MS]|H$', "", cigars))) 
    m_qry[!is_plus] <- shift(m_qry[!is_plus], hclipr[!is_plus])  
  } else {
    m_qry[!is_plus] <- shift(m_qry[!is_plus], hclipl[!is_plus])
  }
  
  plus_min <- c(which(is_plus), which(!is_plus))
  plus_min <- rep(plus_min, lapply(m_qry, length))
  ord <- order(plus_min)

  xs <- mapply(seq, unlist(start(m_qry[is_plus])), unlist(end(m_qry[is_plus])), 
               SIMPLIFY = FALSE)
  if (two_strands){
    xs <- c(xs, mapply(function(x,y) rev(seq(x,y)), unlist(start(m_qry[!is_plus])), 
          unlist(end(m_qry[!is_plus])), SIMPLIFY = FALSE))[ord]
  } else {
    xs <- c(xs, mapply(function(x,y) seq(x,y), unlist(start(m_qry[!is_plus])), 
          unlist(end(m_qry[!is_plus])), SIMPLIFY = FALSE))[ord]
  }
  
  ys <- mapply(seq, unlist(start(m_ref)), unlist(end(m_ref)), SIMPLIFY = FALSE)[ord]
  chrs <- seqnames(genomic_locs)[plus_min][ord]
  
  # Introduce gaps for aligned segments separated by more than max_gap
  y_lns <- sapply(ys, length)
  y_sums <- cumsum(y_lns)
  n_segs <- length(ys)
  
  offsets <- c(0, unlist(ys)[y_sums[-n_segs] +1] - (unlist(ys)[y_sums[-n_segs]]+1))
  offsets[offsets > max_gap] <- max_gap # clip large gaps
  chr_chgs <- cumsum(chrs@lengths)+1 
  offsets[chr_chgs[chr_chgs <= n_segs]] <- max_gap
  
  # Setup data and plot
  pt_coords <- data.frame(x = unlist(xs), ylabs = unlist(ys))
  pt_coords$y = 1:nrow(pt_coords) + rep(cumsum(offsets), y_lns)
  
  ymax <- pt_coords$y[y_sums]
  ymin <- na.omit(c(1, pt_coords$y[y_sums +1]))
  
  # Make boxes around single cigar segments so indels can be distinguished from chimeras
  m_qry_range <- range(m_qry)
  m_qry_start <- unlist(start(m_qry_range))
  m_qry_end <- unlist(end(m_qry_range))
  
  box_ymins <- pt_coords$y[sapply(m_qry_start, function(x) max(which(pt_coords$x == x)))]
  box_ymaxs <- pt_coords$y[sapply(m_qry_end, function(x) min(which(pt_coords$x == x)))]
  
  box_coords <- data.frame(xmin = m_qry_start, xmax = m_qry_end, chrs = seqnames(genomic_locs),
                           #ymin = pt_coords[pt_coords$x %in% m_qry_start, "y"],
                           #ymax = pt_coords[pt_coords$x %in% m_qry_end, "y"])
                           ymin  = box_ymins, ymax = box_ymaxs)
                           
  gap_starts <- y_sums[which(offsets >= max_gap) -1]

  # chr_box_coords are coordinates for shading to indicate gaps
  chr_box_coords <- data.frame(ymin = pt_coords[gap_starts, "y"], 
                               ymax = pt_coords[gap_starts + 1, "y"])
   
  zip <- c(1, c(t(matrix(c(gap_starts, gap_starts + 1), ncol = 2))), nrow(pt_coords))
  
  tick_pos <- unlist(lapply(seq(2, length(zip), by = 2), function(i){
                     diff <- pt_coords[zip[i],"y"]-pt_coords[zip[i-1],"y"] + 1
                     unique(c(seq(0,diff, by = tick_sep), diff)
                        + pt_coords[zip[i-1], "y"])}))
  
  tick_labs <- unlist(lapply(seq(2, length(zip), by = 2), function(i){
                     diff <- pt_coords[zip[i],"y"]-pt_coords[zip[i-1],"y"] + 1
                      unique(c(seq(0,diff, by = tick_sep), diff)
                        + pt_coords[zip[i-1], "ylabs"])}))
  
  xbreaks <- seq(min(pt_coords$x), max(pt_coords$x), by = tick_sep)
  
  # Plotting
  p <- ggplot(pt_coords, aes(x=x,y=y)) + geom_point() +
       geom_rect(data = box_coords, aes(xmin = xmin, xmax = xmax, 
                  ymin = ymin, ymax = ymax, x = NULL, y = NULL, fill = chrs, colour = chrs), 
                  alpha = 0.25) + 
       scale_y_continuous(expand = c(0,0), breaks = tick_pos, labels = tick_labs)+
       xlab("Read location") + ylab("Chromosomal location") +
       guides(fill = guide_legend(title = legend_title), 
              colour = guide_legend(title = legend_title)) +
       theme_bw() + theme(axis.title.y=element_text(vjust = 2, size = title_size), 
                          axis.title.x=element_text(vjust = -1, size = title_size),
                          axis.text.y=element_text(size=text_size),
                          axis.text.x=element_text(size=text_size),
                          plot.margin = unit(c(1, 1, 1, 1), "lines")) 
  if (nrow(chr_box_coords) > 0){                         
    p <- p + geom_rect(data = chr_box_coords, aes(ymin = ymin, 
                       ymax = ymax, y = NULL, x = NULL), 
                       xmin = min(pt_coords$x), xmax = max(pt_coords$x), 
                       fill = "gray", alpha = 0.2, colour = "gray", linetype = "dotted") +
         scale_x_continuous(expand = c(0,0), breaks = xbreaks) 
  }
  p
}