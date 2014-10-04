library(GenomicAlignments)
library(GenomicRanges)
library(ggplot2)
library(gridExtra)


plotChimeras <- function(chimera_alns, max_gap = 10, tick_sep = 20, text_size = 12,  
                         title_size = 16, legend_title = "Chromosome"){
  # max_gap: if genomic segments are separated by more than max_gap, 
  #          a gap in the y-axis will be introduced
  
  cigars <- cigar(chimera_alns)
  genomic_locs <- as(chimera_alns, "GRanges")
  
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
  hclipr <- as.numeric(gsub("^$", 0, gsub('.*[MS]|H$', "", cigars))) 
  is_plus <- as.vector(strand(genomic_locs) == "+")
  m_qry_plus <- shift(m_qry[is_plus], hclipl[is_plus])
  m_qry_minus <- shift(m_qry[!is_plus], hclipr[!is_plus])
  plus_min <- c(which(is_plus), which(!is_plus))
  ord <- order(plus_min)

  xs <- mapply(seq, start(m_qry_plus), end(m_qry_plus))
  xs <- c(xs, mapply(function(x,y) rev(seq(x,y)), start(m_qry_minus), end(m_qry_minus),
          SIMPLIFY = FALSE))[ord]
  ys <- mapply(seq, start(m_ref[plus_min]), end(m_ref[plus_min]))[ord]
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
  box_coords <- data.frame(xmin = sapply(xs, min), xmax = sapply(xs, max), ymax = ymax,
                           ymin = ymin, chrs = as.factor(chrs))
   
  gap_starts <- y_sums[which(offsets >= max_gap) -1]
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
                          
  p <- p + geom_rect(data = chr_box_coords, aes(xmin = min(pt_coords$x), 
                    xmax = max(pt_coords$x), ymin = ymin, ymax = ymax, y = NULL, x = NULL), 
                    fill = "gray", alpha = 0.2, colour = "gray", linetype = "dotted") +
       scale_x_continuous(expand = c(0,0), breaks = xbreaks) 
  p 
  
  
}