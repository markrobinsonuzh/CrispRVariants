#'@title Plot alignments, frequencies and location of target sequence
#'@rdname plotVariants
#'@export
setGeneric("plotVariants", function(obj, ...) {
  standardGeneric("plotVariants")})

#'@rdname plotVariants
#'@param obj 
#'@param txdb GenomicFeatures:TxDb object
#'@param plotAlignments.args Extra arguments for plotAlignments
#'@param plotFreqHeatmap.args Extra arguments for plotFreqHeatmap
#'@param add.chr If target chromosome does not start with "chr", e.g. 
#'"chr5", add the "chr" prefix.  (Default:TRUE)
#'@param ... extra arguments for plot layout
setMethod("plotVariants", signature("CrisprSet"),  
          function(obj, ..., txdb, add.chr = TRUE, 
                   plotAlignments.args = list(),
                   plotFreqHeatmap.args = list()){
            
  if(!(class(txdb) == "TxDb" | class(txdb) == "TranscriptDb") ){
    stop("txdb should be a (GenomicFeatures) transcript database object")
  }
  target <- obj$target
  if (add.chr == TRUE){
    # If adding "chr" to target chromosomes matches txdb chromosomes, do so
    target_levels <- seqlevels(target)
    txdb_levels <- seqlevels(txdb)
    wchr <- paste0("chr", target_levels)
    idxs <- wchr %in% txdb_levels
    target_levels[idxs] <- wchr[idxs]  
    target <- renameSeqlevels(target, target_levels)
  }
  gene_p <- annotateGenePlot(txdb, target)
  plotAlignments.args$obj = obj
  aln_p <- do.call(plotAlignments, plotAlignments.args)
  plotFreqHeatmap.args$obj = obj
  heat_p <- do.call(plotFreqHeatmap, plotFreqHeatmap.args) 
  result <- arrangePlots(gene_p, aln_p, heat_p, ...)    
  return(result)
})

#'@title Arrange plots for plotVariants:CrisprSet
#'@description Arranges 3 plots in two rows.  The vertical margins of the
#'left.plot and right.plot constrained to be equal 
#'@param top.plot ggplot grob, placed on top of the figure, spanning the figure width
#'@param left.plot ggplot, placed in the second row on the left
#'@param right.plot ggplot, placed in the second row on the right.  
#'y-axis labels are removed.
#'@param fig.height Actual height for the figure. If not provided,
#'figure height is the sum of the row.ht.ratio (Default: NULL)
#'@param row.ht.ratio Vector specifying row height ratio (Default: c(1,6))
#'@param col.wdth.ratio  Vector specifying column width ratio (Default: c(2, 1))
#'@param left.plot.margin Unit object specifying margins of left.plot.  
#'Margins of right.plot are constrained by the left.plot.
arrangePlots <- function(top.plot, left.plot, right.plot, fig.height = NULL,
                         col.wdth.ratio  = c(2, 1), row.ht.ratio = c(1,6), 
                         left.plot.margin = unit(c(0.25,0,10,0.5), "lines")){
        
  # Set the size ratio of the top and bottom rows
  plot_hts <- if (is.null(fig.height)){ row.ht.ratio 
  }else { fig.height/sum(row.ht.ratio)*row.ht.ratio } 
  
  # Remove y-axis labels from right plot
  right.plot <- right.plot + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank())            
  
  # Adjust margins of left.plot
  left.plot <- left.plot + theme(plot.margin = left.plot.margin)            
 
  # Convert plots to grobs, lock plot heights
  p2 <- ggplotGrob(left.plot)  
  p3 <- ggplotGrob(right.plot)
  p3$heights <- p2$heights  
    
  # Return arranged plots
  return(grid.arrange(top.plot, arrangeGrob(p2, p3, ncol = 2, widths = col.wdth.ratio), 
                      nrow = 2, heights = plot_hts, newpage = FALSE))
}

#'@title Plots gene structure and annotates with target location
#'@description Uses ggbio to plot the gene structure, annotates this with the
#'target location
#'@param txdb A GenomicFeatures:TxDb object
#'@param target Location of target (GRanges)
#'@param target.colour Colour of box indicating targt region
#'@param target.size Thickness of box indicating target region
#'@param gene.text.size Size for figure label
#'@param panel.margin Unit object, margin size
annotateGenePlot <- function(txdb, target, target.colour = "red", target.size = 1, 
                             gene.text.size = 20, 
                             panel.margin = unit(c(0.25,0.25,0.25,0.25), "lines")){
  # Make the gene plot
  stopifnot(require(ggbio))
  genes <- genes(txdb)
  wh <- genes[findOverlaps(genes, target, ignore.strand = TRUE)@queryHits]  
  
  if (length(wh) == 0){
    p1 <- grob()
  } else{
    gene_id <- mcols(wh)$gene_id
    cat("Creating transcript plot with ggbio\n")
    p1 <- ggbio::autoplot(txdb, wh, label = FALSE)
    
    #Pull off the y limits from the transcript plot
    yranges <- ggplot_build(p1)$panel$ranges[[1]]$y.range 
    xranges <- ggplot_build(p1)$panel$ranges[[1]]$x.range
    target_df <- data.frame(xmin = start(target), xmax = end(target), 
                            ymin = yranges[1], ymax = yranges[2])
    
    p1 <-  p1 + geom_rect(data = target_df, 
                          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                          colour = target.colour, fill = NA, size = target.size) + 
      theme_minimal() + ggtitle(gene_id) + 
      theme(axis.text.x = element_text(size = gene.text.size), 
            text = element_text(size = gene.text.size),
            panel.margin = panel.margin) 
    
    p1 <- as.list(attributes(p1))$ggplot
    p1 <- ggplotGrob(p1)    
  }
  return(p1)
}