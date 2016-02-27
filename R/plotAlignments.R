#'@title Plot alignments with respect to a reference sequence
#'@rdname plotAlignments
#'@param obj The object to be plotted
#'@param ... Additional arguments
#'@return A ggplot2 figure
#'@seealso \code{\link{seqsToAln}}, \code{\link[ggplot2]{ggplot}}
#'@author Helen Lindsay
#'@export
setGeneric("plotAlignments", function(obj, ...) {
  standardGeneric("plotAlignments")})

#'@rdname plotAlignments
#'@description (signature("CrisprSet")) Wrapper for CrisprSet$plotVariants.
#'Optionally filters a CrisprSet frequency table, then plots variants.
#'More information in \code{\link[CrispRVariants]{CrisprSet}}
#'@param min.freq i (%) only plot variants with frequency >= i% in at least
#' one sample (default: 0, i.e no frequency cutoff)
#'@param min.count i (integer) only plot variants with count >= i in at least
#' one sample (default: 0, i.e no count cutoff)
#'@param top.n (integer) Plot only the n most frequent variants
#' (default: 50)
#'@param renumbered If TRUE, the x-axis is numbered with respect to the target
#' (default: TRUE)
#'@param add.other Add a blank row labelled "Other" to the plot, for combining
#'with plotFreqHeatmap (default: TRUE (signature "CrisprSet") FALSE (signature "matrix"))
#'@examples
#'#Load a CrisprSet object and plot
#'data("gol_clutch1")
#'plotAlignments(gol)
setMethod("plotAlignments", signature("CrisprSet"),
          function(obj, ..., min.freq = 0, min.count = 1,
                   top.n = 50, renumbered = obj$pars[["renumbered"]],
                   add.other = TRUE) {

            plot_obj <- obj$plotVariants(min.freq = min.freq,
                                         min.count = min.count,
                                         top.n = top.n,
                                         renumbered = renumbered,
                                         add.other = add.other, ...)

            return(plot_obj)
          })


#'@title Plots pairwise alignments
#'@description (signature("DNAString"))  Plots a set of pairwise alignments to a reference sequence.
#'Alignments should all be the same length as the reference sequences.
#'This is achieved by removing insertions with respect to the reference,
#'see \code{\link[CrispRVariants]{seqsToAln}}.
#'Insertions are indicated by symbols in the plot and a table showing the
#'inserted sequences below the plot.  The default options are intended for a
#'figure 6-8 inches wide, with figure height best chosen according to the number
#'of different variants and insertions to be displayed.
#'
#'@param alns A named character vector of aligned sequences, with insertions removed
#'@param ins.sites A table of insertion_sites, which must include cols
#'named "start", "cigar" and "seq", for the start of the insertion in the
#'corresponding sequence
#'@param highlight.pam should location of PAM with respect to the target site be
#'indicated by a box? (Default: TRUE)  If TRUE, and pam.start and pam.end are not
#'supplied, PAM is inferred from target.loc
#'@param show.plot  Should the plot be displayed (TRUE) or just returned as
#'a ggplot object (FALSE).  (Default: FALSE)
#'@param target.loc The location of the zero point / cleavage location.  Base n, where
#'the zero point is between bases n and n+1
#'@param pam.start  The first location of the PAM with respect to the reference.
#'@param pam.end    The last location of the PAM with respect to the reference.
#'Default is two bases after the pam.start
#'@param ins.size   The size of the symbols representing insertions within the plot.
#'@param legend.cols  The number of columns in the legend.  (Default:3)
#'@param xlab   A title for the x-axis (Default: NULL)
#'@param xtick.labs Labels for the x-axis ticks (Default: NULL)
#'@param xtick.breaks Locations for x-axis tick breaks (Default: NULL)
#'@param plot.text.size The size of the text inside the plot
#'@param axis.text.size The size of the axis labels
#'@param legend.text.size The size of the legend labels
#'@param highlight.guide  Should the guide be indicated by a box in
#'the reference sequence?  (Default: TRUE)
#'@param guide.loc  The location of the guide region to be highlighted,
#' as an IRanges object. Will be inferred from target.loc if
#' highlight.guide = TRUE and no guide.loc is supplied, assuming the guide
#' plus PAM is 23bp (Default: NULL)
#'@param tile.height  The height of the tiles within the plot. (Default: 0.55)
#'@param max.insertion.size  The maximum length of an insertion to be shown in the
#'legend.  If max.insertion.size = n, an insertion of length m > n will
#'be annotated as "mI" in the figure.  (Default: 20)
#'@param min.insertion.freq  Display inserted sequences with frequency at least x
#'amongst the sequences with an insertion of this size and length (Default: 5)
#'@param line.weight  The line thickness for the vertical line indicating the
#'zero point (cleavage site) and the boxes for the guide and PAM.  (Default: 1)
#'@param legend.symbol.size The size of the symbols indicating insertions
#'in the legend.  (Default: ins.size)
#'@rdname plotAlignments
setMethod("plotAlignments", signature("DNAString"),
  function(obj, ..., alns, ins.sites, highlight.pam = TRUE, show.plot = FALSE,
           target.loc = 17, pam.start = NA, pam.end = NA,
           ins.size = 2, legend.cols = 3, xlab = NULL, xtick.labs = NULL,
           xtick.breaks = NULL, plot.text.size = 2, axis.text.size = 8,
           legend.text.size = 6, highlight.guide=TRUE, guide.loc = NULL,
           tile.height = 0.55, max.insertion.size = 20, min.insertion.freq = 5,
           line.weight = 1, legend.symbol.size = ins.size, add.other = FALSE){

  # Insertion locations are determined by matching ins.sites$cigar with names(alns)
  ref <- obj
  m <- transformAlnsToLong(ref, alns, add.other = add.other)
  m <- setDNATileColours(m)
  nms <- m$Var1[1:(length(alns) + 1)]

  p <- makeAlignmentTilePlot(m, ref = ref, xlab = xlab,
          plot.text.size = plot.text.size, axis.text.size = axis.text.size,
          xtick.labs = xtick.labs, xtick.breaks = xtick.breaks,
          tile.height = tile.height)

  # Colours and shapes for the insertion markers and tiles
  shps <- c(21,23,25)
  clrs <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00",
            "#CC79A7","#332288","#88CCEE","#44AA99","#117733","#999933",
            "#DDCC77","#661100","#CC6677","#882255", "#AA4499")

  # Add line for the cut site
  p <- p + geom_vline(xintercept = target.loc + 0.5,
                      colour = "black", size = line.weight)

  top_row <- ifelse(add.other, length(nms) + 1, length(nms))

  if (highlight.guide){
    ymin = top_row - (tile.height / 2 + 0.25)
    ymax = top_row + (tile.height / 2 + 0.25)

    if (is.null(guide.loc)){
      xmin <- target.loc - 16.5
      xmax <- xmin + 23
    } else {
      xmin <- start(guide.loc) - 0.5
      xmax <- end(guide.loc) + 0.5
    }

    guide_df <- data.frame(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
    p <- p + geom_rect(data = guide_df,
                       aes_q(xmin = quote(xmin), xmax = quote(xmax),
                             ymin = quote(ymin), ymax = quote(ymax),
                             color = "black"),
                       size = line.weight, fill = "transparent")
  }

  # If pam_loc is given, highlight the pam in the reference
  #  - 0.5 for tile boundaries not centres
  if (isTRUE(highlight.pam)){
    if (! is.na(pam.start)){
      if (is.na(pam.end)){
        pam.end <- pam.start + 2
      }
    } else {
      # Infer from the target location
      pam.start <- target.loc + 4
      pam.end <- target.loc + 6
    }

    pam_df <- data.frame(xmin = pam.start - 0.5, xmax = pam.end + 0.5,
                         ymin = top_row - (tile.height / 2 + 0.2),
                         ymax = top_row + (tile.height / 2 + 0.2))
    p <- p + geom_rect(data=pam_df,
                       aes_q(xmin = quote(xmin), xmax = quote(xmax),
                             ymin = quote(ymin), ymax = quote(ymax)),
                       color = "black", size = line.weight, fill = "transparent")
    #p <- p + annotation_custom(grob = textGrob("PAM", gp = gpar(cex = 3)),
    #                           xmin = pam_df$xmin, xmax = pam_df$xmax,
    #                           ymin = pam_df$ymin + 1, ymax = pam_df$ymax + 1)

  }

  # Make a data frame of insertion locations
  ins_ord <- match(ins.sites$cigar, nms)

  if (nrow(ins.sites) > 0 & length(na.omit(ins_ord)) > 0){
    ins_points <- data.frame(x = ins.sites[!is.na(ins_ord),"start"] - 0.5,
                             y = na.omit(ins_ord) + 0.45,
                             seq = ins.sites[!is.na(ins_ord),"seq"],
                    count = as.integer(ins.sites[!is.na(ins_ord),"count"]))

    ins_points <- aggregate(ins_points$count,
                            by = as.list(ins_points[,c(1:3)]), FUN = sum)
    colnames(ins_points) <- c("x","y","seq", "count")

    ## Merge multiple insertions at single plotting location, format to fixed width
    xy_locs <- paste(ins_points$x, ins_points$y, sep = "_")

    # Remove low frequency alleles
    temp <- rowsum(ins_points$count, xy_locs)
    totals <- as.vector(temp)
    names(totals) <- rownames(temp)
  
    temp <- ! duplicated(xy_locs)
    wdths <- nchar(as.character(ins_points[temp, "seq"]))
    
    omit <- (ins_points$count/totals[xy_locs])*100 < min.insertion.freq
    ins_points[omit,"seq"] <- NA

    names(wdths) <- xy_locs[temp]
    splits <- split(ins_points$seq, xy_locs)
    
    # Collapse sequences longer than max.insertion.size
    x <- lapply(seq_along(splits), function(i){
      result <- as.character(na.omit(splits[[i]]))
      wdth <- wdths[[names(splits)[i]]]
      
      if (length(result) != 1){
        # If there are multiple sequences
        if (wdth > max.insertion.size | length(result) == 0){
          result <- sprintf("%sI (%s common alleles)", wdth, length(result))
        } else {  
          # Collapse, one sequence per line
          result <- paste(result, collapse = ",\n")
        }
      } else {
          # If only one sequence
          if (wdth > max.insertion.size) result <- sprintf("%sI", wdth)
      }
      result
    })
    names(x) <- names(splits)

    new_seqs <- unlist(x)[unique(xy_locs)]
    max_seq_ln <- max(sapply(gsub("\n.*", "", new_seqs), nchar)) + 3
    new_seqs <- sprintf(paste0("%-",max_seq_ln,"s"), new_seqs)
    
    ins_points <- ins_points[!duplicated(ins_points[,c("x","y")]),]
    ins_points$seq <- new_seqs
    
    # Specify colours and shapes for insertion symbols
    ins_points$shapes <- as.factor(1:nrow(ins_points))
    ins_points$colours <- as.factor(rep(clrs, 3)[1:nrow(ins_points)])
    fill_clrs <- rep(clrs, 3)[1:nrow(ins_points)]
    fill_shps <- rep(shps, 17)[1:nrow(ins_points)]

    legend_nrow <- ceiling(nrow(ins_points)/ legend.cols)

    # Indicate insertions
    p <- p + geom_point(data = ins_points,
                        aes_q(x = quote(x), y = quote(y),
                              shape = quote(shapes), fill = quote(colours)),
                        colour = "#000000", size = ins.size)  +
      scale_fill_identity() +
      scale_shape_manual(name = "", values = fill_shps, breaks = ins_points$shapes,
                         labels = ins_points$seq) +
      guides(shape = guide_legend(nrow = legend_nrow,
                            override.aes = list(fill = fill_clrs,
                                                size = legend.symbol.size)))
    p <- p + theme(legend.key = element_blank(),
                   legend.text = element_text(size = legend.text.size),
                   legend.margin=grid::unit(0.2,"cm"))

  } else{
    p <- p + scale_fill_identity()
  }

  if (isTRUE(show.plot)){
    print(p)
  }
  return(p)
})


#'@title Transform data for plotting
#'@description Orders and transforms a reference sequence and a set of aligned sequences
#'into long format, i.e. one observation (tile position) per row.  Used internally by
#'\code{\link[CrispRVariants]{plotAlignments}}.
#'@param ref The reference sequence
#'@param alns Character vector of aligned sequences
#'@param add.other Add a blank row labelled "Other" (Default: FALSE)
#'@return A matrix of characters and plotting locations
#'@author Helen Lindsay
transformAlnsToLong <- function(ref, alns, add.other = FALSE){
  # Reverse alignment order, as ggplot geom_tile plots bottom up
  aln_chrs <- strsplit(c(rev(alns), Reference = as.character(ref)), "")

  # Add a blank row which will get white tiles
  if (isTRUE(add.other)){
    aln_chrs <- c(list(rep("",length(aln_chrs[[1]]))), aln_chrs)
    names(aln_chrs)[[1]] <- "Other"
  }

  # Test that all alignments and reference have the same length
  if (! length(unique(lapply(aln_chrs, length))) == 1){
    stop("The reference sequence and the alignments must all have the same length")
  }

  temp <- t(as.data.frame(aln_chrs))
  rownames(temp) <- names(aln_chrs)
  m <- reshape2::melt(temp)

  m
}


#'@title Sets colours for plotting aligned DNA sequences.
#'@description Sets tile colours for \code{\link[CrispRVariants]{plotAlignments}} with a
#'DNA alphabet.  Colour names must be valid.
#'@param m A matrix with a column named "value" of the characters at each tile position.
#'@return A matrix with additional columns specifying tile and text colours
#'@author Helen Lindsay
setDNATileColours <- function(m){
  ambig_codes <- c('K','M','R','Y','S','W','B','V','H','D')
  ambig <- which(! m$value %in% c("A", "C", "T", "G", "N", "-"))

  m$value <- as.character(m$value)
  m$value <- factor(m$value, levels = c(c("A", "C", "T", "G", "N", "-", ""),
                                        ambig_codes))
  m$isref <- as.character(ifelse(m$Var1 == "Reference", 1, 0.75))
  m_cols <- c(c("#4daf4a", "#377eb8", "#e41a1c", "#000000", "#CCCCCC","#FFFFFF",
                "#FFFFFF"), rep("#CCCCCC", length(ambig_codes)))
  names(m_cols) <- c(c("A", "C", "T", "G", "N","-"), ambig_codes)
  m$cols <- m_cols[m$value]
  m$text_cols <- ifelse(m$cols == "#000000" & m$isref == 1, "#FFFFFF", "#000000")
  return(m)
}


makeAlignmentTilePlot <- function(m, ref, xlab, plot.text.size, axis.text.size,
                                  xtick.labs, xtick.breaks, tile.height){
  alpha_v <- c(0.5,1)
  # Change alpha values if only plotting the reference
  if (length(unique(m$Var1)) == 1) alpha_v <- 1
  # Plot aligned sequences
  p <- ggplot(m) +
    geom_tile(aes_q(x = quote(Var2), y = quote(Var1),
                    fill = quote(cols), alpha = quote(isref)),
                    height = tile.height) +
    geom_text(aes_q(x = quote(Var2), y = quote(Var1),
                    fill = quote(cols), label = quote(value),
                    colour = quote(text_cols)), size = plot.text.size) +
    scale_alpha_manual(values = alpha_v, guide = "none") +
    ylab(NULL) + xlab(xlab) + scale_colour_identity() +
    theme_bw() + theme(axis.text.y = element_text(size = axis.text.size),
                       axis.text.x = element_text(size = axis.text.size),
                       axis.title.x = element_text(vjust = -0.5),
                       legend.position = "bottom")

  if (is.null(xtick.labs)){
    # expand is the distance from the axis, multiplicative + additive
    p <- p + scale_x_continuous(expand = c(0,0.25))
  } else {
    if (is.null(xtick.breaks)) {
      stopifnot(length(xtick.labs == nchar(ref)))
      p <- p + scale_x_continuous(expand = c(0,0.25), breaks = 1:nchar(ref),
                                  labels = xtick.labs)
    } else {
      p <- p + scale_x_continuous(expand = c(0,0.25), breaks = xtick.breaks,
                                  labels = xtick.labs)
    }

  }
  return(p)
}
