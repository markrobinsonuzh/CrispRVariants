#'@title Plot a table of counts with colours indicating frequency
#'@rdname plotFreqHeatmap
#'@export
setGeneric("plotFreqHeatmap", function(obj, ...) {
  standardGeneric("plotFreqHeatmap")})


#'@rdname plotFreqHeatmap
#'@param obj A matrix of counts with rows = feature, columns = sample
#'@param col.sums Include a row of column totals at the top of the
#'plot (Default: TRUE)
#'@param header Alternative column titles, e.g. column sums
#'for the unfiltered data set (Default: NULL)
#'@param group Grouping factor for columns.  If supplied, columns are
#'ordered to match the levels  (Default: NULL)
#'@param group.colours Colours for column groups, should match levels of "group".
#'If "NULL", groups are coloured differently (Default: NULL)
#'@param as.percent  Should colours represent the percentage of reads 
#'per sample (TRUE) or the actual counts (FALSE)?  (Default: TRUE)   
#'@param x.axis.title A title for the x-axis.  (Default: NULL)
#'@param x.size Font size for x-labels (Default: 16)
#'@param y.size Font size for y-labels (Default: 16)
#'@param x.angle  Angle for x-labels (Default: 90, i.e. vertical)
#'@param legend.text.size Font size for legend (Default: 16)
#'@param plot.text.size Font size counts within plot (Default: 8)
#'@param line.width Line thickness of title box
#'@param legend.position The position of the legend (Default: right)
#'@param legend.key.height The height of the legend key, as a "unit" object.
#'(See \code{\link[grid]{unit}}).
#'@param ... additional arguments
#'@return The ggplot2 plot of the variant frequencies 
setMethod("plotFreqHeatmap", signature("matrix"),  
          function(obj, ..., col.sums = TRUE, header = NULL, group = NULL,
                   group.colours = NULL, as.percent = TRUE, x.axis.title = NULL,
                   x.size = 6, y.size = 8, x.angle = 90, legend.text.size = 6,
                   plot.text.size = 2, line.width = 1, legend.position = "right",
                   legend.key.height = grid::unit(2, "lines")) {            
  
  # Make space for totals to be added
  if (length(header) == ncol(obj)){
    col.sums <- TRUE
  } else if (length(header) > 0){
    stop("Header length should equal to the number of columns of obj")
  }
  
  if (col.sums == TRUE){
    obj <- rbind(Total = rep(NA, ncol(obj)), obj)
  }

  # If a sample group is supplied, reorder the columns of counts  
  if (! is.null(group)){
    if (! class(group) == "factor") group <- factor(group, levels = unique(group)) 
    
    obj <- obj[,order(group), drop = FALSE] 
    group <- group[order(group)]
    
    if (is.null(group.colours)){
      clrs <- c("#332288","#661100","#117733","#882255","#D55E00", 
                "#0072B2","#AA4499","#009E73","#56B4E9","#CC79A7",
                "#44AA99","#999933","#CC6677","#E69F00","#88CCEE")
      clrs <- clrs[group]
    } else { 
      clrs <- group.colours[group] 
    }   
  }    
  
  counts <- reshape2::melt(obj)  
  colnames(counts) <- c("Feature", "Sample","Count")
  counts$Feature <- factor(counts$Feature, levels = rev(levels(counts$Feature)))
  
  # Create coloured tile background
  if (as.percent == TRUE){
    m <- apply(obj, 2, function(x) x/sum(na.omit(x)))  
    m <- melt(m)
    colnames(m) <- c("Feature", "Sample","Percentage")  
    m$Feature <- factor(m$Feature, levels = rev(levels(m$Feature)))
    g <- ggplot(m, aes(x = Sample, y = Feature, fill = Percentage)) + geom_tile()
  }
  else{    
    g <- ggplot(counts, aes(x = Sample, y = Feature, fill = Count)) + geom_tile()
  }
  
  # Add the count numbers to the boxes,
  # if totals are included, add bold boxes to indicate these
  counts$ff <- "plain"
  box_coords <- data.frame(xmin = integer(), xmax=integer(),
                           ymin=integer(), ymax=integer())
  box_row <- 1
  xranges <- ggplot_build(g)$panel$ranges[[1]]$x.range  
  yranges <- ggplot_build(g)$panel$ranges[[1]]$y.range  
  
  if (col.sums == TRUE){ 
    idxs <- which(is.na(counts$Count))
    if (length(header) != ncol(obj)){
      header <- colSums(obj, na.rm = TRUE)
    }
    counts$Count[idxs] <- header
    counts$ff[idxs] <- "bold"
    box_coords[box_row,] <- c(min(xranges), max(xranges), 
                             nrow(obj) -0.5,max(yranges))  
  }
  
  if (nrow(box_coords) > 0){
    g <- g + geom_rect(data=box_coords, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax,
                                            ymax = ymax, x = NULL, y = NULL),
                       color = "black", size = line.width, fill = "transparent")   
  }
  
  # Plot the counts
  g <- g + geom_text(data = counts, aes(label = Count, fill = NULL, fontface = ff), 
                     size = plot.text.size)
  
  # Colour the boxes - white for 0, darkred for highest
  hmcols<-colorRampPalette(c("white","gold","orange","orangered","red", "darkred"))(50) 
  
  # No colour bar legend if the max count is 1
  if (max(counts$Count) > 1){
    g <- g + scale_fill_gradientn(colours = hmcols, na.value = "white", guide = "legend") 
  } else {
    g <- g + scale_fill_gradientn(colours = hmcols, na.value = "white", guide = "none")
  }
  
  # Set plot labels
  g <- g + ylab(NULL) + xlab(x.axis.title) + theme_bw() +
    theme(axis.text.x = element_text(size = x.size, angle = x.angle, hjust = 1),
          axis.text.y = element_text(size = y.size),
          legend.text = element_text(size = legend.text.size),
          legend.title = element_text(size = legend.text.size),
          legend.position = legend.position,
          legend.key.height = legend.key.height)         
  
  # Colour xlabels by group, make sure not to replace original args
  if (! is.null(group)){
    g <- g + theme(axis.text.x=element_text(colour= clrs, 
                   size = x.size, angle = x.angle, hjust = 1))
  }
  
  return(g)
})


#'@rdname plotFreqHeatmap
#'@param top.n  Show the n top ranked variants.  Note that if the nth and n+1th 
#'variants have equal rank, they will not be shown.   (Default: 50)
#'@param min.freq i (%) only plot variants with frequency >= i% in at least
#' one sample (default: 0, i.e no frequency cutoff)
#'@param min.count i (integer) only plot variants with count >= i in at least
#' one sample (default: 0, i.e no count cutoff)
#'@param type Plot either "counts" or "proportions"
#'@examples
#'#Load a CrisprSet object for plotting
#'data("gol_clutch1")
#'
#'# Plot the frequency heatmap
#'plotFreqHeatmap(gol)
setMethod("plotFreqHeatmap", signature("CrisprSet"),  
          function(obj, ..., top.n = 50, min.freq = 0, min.count = 0, 
                   type = c("counts", "proportions")) {
  
  result <- obj$heatmapCigarFreqs(top.n = top.n, min.freq = min.freq,
                                  min.count = min.count, type = type, ...)          
  return(result)
})      
