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
  
  dots <- list(...)
  annotate_nms = c("target.colour","target.size","gene.text.size", "panel.margin")
  annotate_args = dots[names(dots) %in% annotate_nms]
  dots[annotate_nms] <- NULL
  
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
  annotate_args <- modifyList(list(txdb = txdb, target = target), annotate_args)
  gene_p <- do.call(annotateGenePlot, annotate_args)
  
  plotAlignments.args$obj = obj
  aln_p <- do.call(plotAlignments, plotAlignments.args)
  plotFreqHeatmap.args$obj = obj
  heat_p <- do.call(plotFreqHeatmap, plotFreqHeatmap.args) 
  arrange_args = list(gene_p, aln_p, heat_p)
  result <- do.call(arrangePlots, arrange_args)
  
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
                         left.plot.margin = unit(c(0.1,0,8,0.2), "lines")){
        
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
#'@param plot.title A title for the plot.  If no plot.title is supplied, the title is the 
#'list of gene ids shown (default).  If plot.title == FALSE, the plot will not have a title.
#'@param all.transcripts If  TRUE (default), all transcripts of genes overlapping 
#'the target are shown, including transcripts that do not themselves overlap the target.
#'If FALSE, only the transcripts that overlap the target are shown.
annotateGenePlot <- function(txdb, target, target.colour = "red", target.size = 1, 
                             gene.text.size = 12, 
                             panel.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),
                             plot.title = NULL, all.transcripts = TRUE)){
  
  exByTx <-exonsBy(txdb,"tx")
  utr5 <- fiveUTRsByTranscript(txdb)
  utr3 <- threeUTRsByTranscript(txdb) 
  
  exs <- findOverlaps(guide, exByTx)
  
  if (all.transcripts == TRUE){
    genes <- suppressWarnings(select(txdb, key = unique(genes$GENEID), 
                  keytype = "GENEID", columns = c("GENEID", "TXID", "TXNAME")))
  } else {
    genes <- select(txdb, key = as.character(subjectHits(exs)), 
                    keytype = "TXID", columns = c("GENEID","TXNAME"))
  }
  
  txid <- as.character(genes$TXID)
  # Get UTRs matching each transcript
  all_sections <- lapply(txid, function(tid) {
    u5 <- GRangesList()
    u3 <- GRangesList()
    if (tid %in% names(utr5)){ u5 <- utr5[tid] }
    if (tid %in% names(utr3)){ u3 <- utr3[tid] }
    utrs <- unlist(c(u5, u3))
    split_ex <- disjoin(c(exByTx[[tid]], utrs))
    split_ex$type = "exon"
    split_ex$type[split_ex %in% utrs] <- "utr"
    split_ex
  })

  all_exs <- do.call(c, all_sections)
  min_st <- min(start(all_exs))
  max_end <- max(end(all_exs))
  
  all_exs <- data.frame(start = start(all_exs), 
                        end = end(all_exs), 
                        ts = rep(1:length(all_sections), lapply(all_sections, length)),
                        type = all_exs$type)
  
  colnames(all_exs) <- c("start", "end", "ts", "type")
  
  gene_spans <- do.call(c,unname(lapply(all_sections, range)))
  
  tcks <- unname(quantile(min_st:max_end, seq(1,100, by = 2)*0.01))
  tcks <- c(tcks, start(gene_spans), end(gene_spans))
  tcks <- lapply(gene_spans, function(sp){
    tcks[tcks >= start(sp) & tcks <= end(sp)]  
  })
  
  tcks <- data.frame(tloc = unlist(tcks), ys = rep(1:length(tcks), lapply(tcks, length)))

  all_exs$ymax <- all_exs$ts + 0.3
  all_exs$ymin <- all_exs$ts - 0.3 
  is_utr <- all_exs$type == "utr"
  all_exs$ymax[is_utr] <- all_exs$ts[is_utr] + 0.2
  all_exs$ymin[is_utr] <- all_exs$ts[is_utr] - 0.2 
 
  target_df <- data.frame(xmin = start(target), xmax = end(target), 
                          ymin = 0, ymax = ceiling(max(all_exs$ymax)))
  
  if (is.null(plot.title)){ plot.title <- paste(unique(genes$GENEID), sep = ";")}
  
  p <- ggplot(tcks, aes(x = tloc, y = ys, group = ys)) + geom_line() + 
    geom_point(shape = 62, size = 3) 
  p <- p + geom_rect(data = all_exs, fill = "black", color = "black", 
              aes(x = NULL, y = NULL, group = NULL, 
                  xmin = start, xmax = end, ymin = ymin, ymax=ymax)) 
  p <- p + geom_rect(data = target_df,
                     aes(x = NULL, y = NULL, group = NULL,
                         xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                         colour = target.colour, fill = NA, size = target.size)
  if (! plot.title == FALSE){
    p <- p + ggtitle(plot.title) 
  }
  p <- p + theme_minimal() +
       theme(axis.text.x = element_text(size = gene.text.size), 
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.margin = panel.margin,
          text = element_text(size = gene.text.size),
          axis.ticks.y = element_blank()) +
       ylab(NULL) + xlab(NULL) 
  
  return(p)
}
