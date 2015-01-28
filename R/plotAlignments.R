#'@title Plot alignments with respect to a reference sequence  
#'@rdname plotAlignments
#'@export
setGeneric("plotAlignments", function(obj, ...) {
  standardGeneric("plotAlignments")})

#'@rdname plotAlignments
#'@description Wrapper for CrisprSet$plotVariants.  Optionally filters a 
#'CrisprSet frequency table, then plots variants.  More information in
#'\code{\link[crispRvariants]{CrisprSet}} 
#'@param freq_cutoff i (integer) only plot variants that occur >= i times
#' (default: 0, i.e no frequency cutoff)
#'@param top_n (integer) Plot only the n most frequent variants 
#' (default: plot all)
#'@param renumbered If TRUE, the x-axis is numbered with respect to the target
#' (default: TRUE)
setMethod("plotAlignments", signature("CrisprSet"),  
          function(obj, ..., freq_cutoff = 0, 
                   top_n = nrow(obj$cigar_freqs),
                   renumbered = obj$pars["renumbered"]) {
            plot_obj <- obj$plotVariants(freq_cutoff = freq_cutoff, top_n = top_n, 
                                         renumbered = renumbered, ...)
            
            return(plot_obj)
          })


#'@title Plots pairwise alignments 
#'@description Plots a set of pairwise alignments to a reference sequence.
#'Alignments should all be the same length as the reference sequences.  
#'This is achieved by removing insertions with respect to the reference, 
#'see \code{\link[crispRvariants]{seqsToAln}} for these alignments.
#'Insertions are indicated by symbols in the plot and a table showing the
#'inserted sequences below the plot.
#'
#'@param ref The reference sequence
#'@param alns A named character vector of aligned sequences, with insertions removed
#'@param ins_sites A table of insertion_sites, which must include cols 
#'named "start", "cigar" and "seq", for the start of the insertion in the 
#'corresponding sequence
#'@param highlight_pam should location of PAM with respect to the target site be shown?
#'(Default: TRUE)  If TRUE, and pam_start and pam_end are not supplied, PAM is inferred
#'from target_loc
#'@param show_plot
#'@param target_loc The location of the zero point / cleavage location.  Base n, where 
#'the zero point is between bases n and n+1
#'@param pam_start
#'@param pam_end
#'@param ins_size 
#'@param legend_cols
#'@param xlab
#'@param xtick_labs
#'@param xtick_breaks
#'@param plot_text_size
#'@param axis_text_size
#'@param legend_text_size
#'@param highlight_guide
#'@param guide_loc
#'@param tile_height
#'@return A ggplot figure  
#'@seealso \code{\link[crispRvariants]{seqsToAln}}, \code{\link[ggplot2]}
#'@author Helen Lindsay
#'@rdname plotAlignments
setMethod("plotAlignments", signature("DNAString"),  
  function(obj, ..., alns, ins_sites, highlight_pam = TRUE, show_plot = FALSE, 
           target_loc = 17, pam_start = NA, pam_end = NA, 
           ins_size = 6, legend_cols = 3, xlab = NULL, xtick_labs = NULL,
           xtick_breaks = NULL, plot_text_size = 8, axis_text_size = 16, 
           legend_text_size = 16, highlight_guide=TRUE, guide_loc = NULL,
           tile_height = 0.55){
  
  # WHY PAM LOC, PAM START AND PAM END?????
  # Insertion locations are determined by matching ins_sites$cigar with names(alns)
  ref <- obj
  
  m <- transformAlnsToLong(ref, alns)
  m <- setDNATileColours(m)
  nms <- m$Var1[1:(length(alns) + 1)]

  p <- makeAlignmentTilePlot(m, ref = ref, xlab = xlab, plot_text_size = plot_text_size, 
          axis_text_size = axis_text_size, xtick_labs = xtick_labs,
          xtick_breaks = xtick_breaks, tile_height = tile_height)
  
  # Colours and shapes for the insertion markers and tiles
  shps <- c(21,23,25) 
  clrs <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
            "#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100",
            "#CC6677","#882255", "#AA4499")
    
  # Add line for the cut site
  p <- p + geom_vline(xintercept= target_loc + 0.5, colour = "black", size = 2)# linetype = "dashed",
  
  # Make a data frame of insertion locations 
  ins_ord <- match(ins_sites$cigar, nms)    

  if (nrow(ins_sites) > 0 & length(na.omit(ins_ord)) > 0){ 
    ins_points <- data.frame(x = ins_sites[!is.na(ins_ord),"start"] - 0.5,
                             y = na.omit(ins_ord) + 0.45, 
                             seq = ins_sites[!is.na(ins_ord),"seq"])
    
    # Merge multiple insertions at single plotting location, format to fixed width
    ins_points <- ins_points[!duplicated(ins_points),]
    
    xy_locs <- paste(ins_points$x, ins_points$y, sep = "_")
    seqs <- ins_sites[!is.na(ins_ord),"seq"]
    splits <- split(ins_points$seq, xy_locs)
    x <- lapply(splits, function(x) paste(as.character(x), collapse = ", "))
    new_seqs <- unlist(x)[unique(xy_locs)]
    max_seq_ln <- max(sapply(new_seqs, nchar)) + 3 
    new_seqs <- sprintf(paste0("%-",max_seq_ln,"s"), new_seqs)
    
    ins_points <- ins_points[!duplicated(ins_points[,c("x","y")]),]
    ins_points$seq <- new_seqs
    
    # Specify colours and shapes for insertion symbols
    ins_points$shapes <- as.factor(1:nrow(ins_points))
    ins_points$colours <- as.factor(rep(clrs, 3)[1:nrow(ins_points)])
    fill_clrs <- rep(clrs, 3)[1:nrow(ins_points)]
    fill_shps <- rep(shps, 17)[1:nrow(ins_points)]
    
    legend_nrow <- ceiling(nrow(ins_points)/ legend_cols)
    
    # Indicate insertions   
    p <- p + geom_point(data = ins_points, aes(x = x, y = y, shape = shapes, fill = colours),
                        colour = "#000000", size = ins_size)  +
      scale_fill_identity() + 
      scale_shape_manual(name = "", values = fill_shps, breaks = ins_points$shapes, labels = ins_points$seq) +
      guides(shape = guide_legend(nrow = legend_nrow, override.aes = list(fill = fill_clrs, size = 8)))#,
    p <- p + theme(legend.key = element_blank(), 
                   legend.text = element_text(size = legend_text_size),
                   legend.margin = unit(2, "lines"))
    
  } else{
    p <- p + scale_fill_identity() 
  }
  
  if (highlight_guide){
    ymin = length(nms) - (tile_height / 2 + 0.25)
    ymax = length(nms) + (tile_height / 2 + 0.25)
    
    if (is.null(guide_loc)){
      xmin <- target_loc - 16.5
      xmax <- xmin + 23
    } else {
      xmin <- start(guide_loc) - 0.5
      xmax <- end(guide_loc) + 0.5
    }
    
    guide_df <- data.frame(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
    p <- p + geom_rect(data = guide_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, 
                                            ymax = ymax, color = "black", x = NULL, y = NULL),
                       size = 2, fill = "transparent") 
  }
  
  # If pam_loc is given, highlight the pam in the reference
  #  - 0.5 for tile boundaries not centres
  if (highlight_pam == TRUE){
    if (! is.na(pam_start)){    
      if (is.na(pam_end)){
        pam_end <- pam_start + 2
      }
    } else {
      # Infer from the target location
      pam_start <- target_loc + 4
      pam_end <- target_loc + 6
    }
      
    pam_df <- data.frame(xmin = pam_start - 0.5, xmax = pam_end + 0.5,
                         ymin = length(nms) - (tile_height / 2 + 0.2),
                         ymax = length(nms) + (tile_height / 2 + 0.2))
    p <- p + geom_rect(data=pam_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax,
                                        ymax = ymax, x = NULL, y = NULL),
                       color = "black", size = 3, fill = "transparent") 
    #p <- p + annotation_custom(grob = textGrob("PAM", gp = gpar(cex = 3)), 
    #                           xmin = pam_df$xmin, xmax = pam_df$xmax, 
    #                           ymin = pam_df$ymin + 1, ymax = pam_df$ymax + 1)
    
  }
  if (show_plot == TRUE){
    print(p)
  }
  return(p)
})


#'@title Transform data for plotting
#'@description Orders and transforms a reference sequence and a set of aligned sequences 
#'into long format, i.e. one observation (tile position) per row.  Used internally by 
#'\code{\link[crispRvariants]{plotAlignments}}.
#'@param ref
#'@param alns
#'@author Helen Lindsay
transformAlnsToLong <- function(ref, alns){
  # Reverse alignment order, as ggplot geom_tile plots bottom up  
  aln_chrs <- strsplit(c(rev(alns), Reference = as.character(ref)), "")
  
  # Test that all alignments and reference have the same length
  if (! length(unique(lapply(aln_chrs, length))) == 1){
    stop("The reference sequence and the alignments must all have the same length")
  }
  
  temp <- t(as.data.frame(aln_chrs))
  rownames(temp) <- names(aln_chrs)
  m <- melt(temp)
  return(m)
}

##'@title Format insertions for plotting
##'@description Formats a table of insertions for plotting with \code{\link[crispRvariants]{plotAlignments}}.
##'@author Helen Lindsay
#formatInsertions <- function(){
#}


#'@title Sets colours for plotting aligned DNA sequences.
#'@description Sets tile colours for \code{\link[crispRvariants]{plotAlignments}} with a 
#'DNA alphabet.  Colour names must be valid.  
#'@param m A matrix with a column named "value" of the characters at each tile position. 
#'@author Helen Lindsay
setDNATileColours <- function(m){
  ambig_codes <- c('K','M','R','Y','S','W','B','V','H','D')
  ambig <- which(! m$value %in% c("A", "C", "T", "G", "N", "-"))
  
  m$value <- as.character(m$value)  
  m$value <- factor(m$value, levels = c(c("A", "C", "T", "G", "N", "-"), ambig_codes))  
  m$isref <- as.character(ifelse(m$Var1 == "Reference", 1, 0.75))
  m_cols <- c(c("#4daf4a", "#377eb8", "#e41a1c", "#000000", "#CCCCCC","#FFFFFF"),
              rep("#CCCCCC", length(ambig_codes)))
  names(m_cols) <- c(c("A", "C", "T", "G", "N","-"), ambig_codes)
  m$cols <- m_cols[m$value]
  m$text_cols <- ifelse(m$cols == "#000000" & m$isref == 1, "#FFFFFF", "#000000")
  return(m)
}

#'@title Makes a basic tile plot of aligned sequences
makeAlignmentTilePlot <- function(m, ref, xlab, plot_text_size, axis_text_size,
                                  xtick_labs, xtick_breaks, tile_height){

  # Plot aligned sequences  
  p <- ggplot(m, aes(x = Var2, y = Var1, fill = cols)) +
    geom_tile(aes(alpha = isref), height = tile_height) + 
    geom_text(aes(label = value, colour = text_cols), size = plot_text_size) + 
    scale_alpha_manual(values = c(0.5,1), guide = "none") + 
    ylab(NULL) + xlab(xlab) + scale_colour_identity() + 
    theme_bw() + theme(axis.text.y = element_text(size = axis_text_size), 
                       axis.text.x = element_text(size = axis_text_size), 
                       axis.title.x = element_text(vjust = -0.5), 
                       legend.position = "bottom")
  
  if (is.null(xtick_labs)){   
    # expand is the distance from the axis, multiplicative + additive
    p <- p + scale_x_continuous(expand = c(0,0.25))  
  } else {
    if (is.null(xtick_breaks)) {
      stopifnot(length(xtick_labs == nchar(ref)))
      p <- p + scale_x_continuous(expand = c(0,0.25), breaks = 1:nchar(ref),
                                  labels = xtick_labs) 
    } else {
      p <- p + scale_x_continuous(expand = c(0,0.25), breaks = xtick_breaks,
                                  labels = xtick_labs) 
    }
  }
  return(p)
}
