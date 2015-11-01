#'@title Plot alignments, frequencies and location of target sequence
#'@rdname plotVariants
#'@param obj The object to be plotted
#'@return A ggplot2 plot of the variants
#'@export
setGeneric("plotVariants", function(obj, ...) {
  standardGeneric("plotVariants")})

#'@rdname plotVariants
#'@param txdb GenomicFeatures:TxDb object (default: NULL)
#'@param plotAlignments.args Extra arguments for plotAlignments
#'@param plotFreqHeatmap.args Extra arguments for plotFreqHeatmap
#'@param add.chr If target chromosome does not start with "chr", e.g.
#'"chr5", add the "chr" prefix.  (Default:TRUE)
#'@param ... extra arguments for plot layout
#'@seealso \code{\link{arrangePlots}} for general layout options
#' and \code{\link{annotateGenePlot}} for options relating
#' to the transcript plot.
#'@examples
#'#Load a CrisprSet object for plotting
#'data("gol_clutch1")
#'
#'#Load the transcript db.  This is a subset of the Ensembl Danio Rerio v73 gtf
#'# for the region 18:4640000-4650000 which includes the targeted gol gene
#'
#'library(GenomicFeatures)
#'fn <- system.file("extdata", "Danio_rerio.Zv9.73.gol.sqlite",
#'                  package = "CrispRVariants")
#'txdb <- loadDb(fn)
#'
#'# Plot the variants
#'p <- plotVariants(gol, txdb = txdb)
#'
#'#In the above plot, the bottom margin is too large, the legend is
#'#cut off, and the text within the plots should be larger.
#'#These issues can be fixed with some adjustments:
#'p <- plotVariants(gol, txdb = txdb,
#'                  plotAlignments.args = list(plot.text.size = 4, legend.cols = 2),
#'                  plotFreqHeatmap.args = list(plot.text.size = 4),
#'                  left.plot.margin = grid::unit(c(0.1,0,0.5,1), "lines"))
#'
setMethod("plotVariants", signature("CrisprSet"),
          function(obj, ..., txdb = NULL, add.chr = TRUE,
                   plotAlignments.args = list(),
                   plotFreqHeatmap.args = list()){

  include_txs <- TRUE
  if(!(class(txdb) == "TxDb" | class(txdb) == "TranscriptDb") ){
    if (is.null(txdb)){
      include_txs <- FALSE
    } else{
      stop("txdb should be a (GenomicFeatures) transcript database object")
    }
  }

  dots <- list(...)
  annotate_nms <- c("target.colour", "target.size", "gene.text.size", "panel.margin")
  annotate_args <- dots[names(dots) %in% annotate_nms]
  dots[annotate_nms] <- NULL

  arrange_nms <-  c("fig.height","col.wdth.ratio", "row.ht.ratio", "left.plot.margin")
  arrange_args <- dots[names(dots) %in% arrange_nms]

  target <- obj$target
  if (isTRUE(add.chr) & isTRUE(include_txs)){
    # If adding "chr" to target chromosomes matches txdb chromosomes, do so
    target_levels <- GenomeInfoDb::seqlevels(target)
    txdb_levels <- GenomeInfoDb::seqlevels(txdb)
    wchr <- paste0("chr", target_levels)
    idxs <- wchr %in% txdb_levels
    target_levels[idxs] <- wchr[idxs]
    target <- GenomeInfoDb::renameSeqlevels(target, target_levels)
  }

  if (isTRUE(include_txs)){
    annotate_args <- modifyList(list(txdb = txdb, target = target), annotate_args)
    gene_p <- do.call(annotateGenePlot, annotate_args)
  } else {
    arrange_args[["row.ht.ratio"]] <- c(0,1)
    gene_p <- grid::grid.rect(gp=grid::gpar(col="white"), draw = FALSE)
    no_ins <- nrow(obj$insertion_sites) == 0
    if (no_ins & ! "left.plot.margin" %in% names(arrange_args)){
      # Sample names tend to clip if there are no insertions, so increase default
      args[["left.plot.margin"]] <- grid::unit(c(0.1,0,3,0.5), "lines")
    }
  }

  plotAlignments.args$obj = obj
  aln_p <- do.call(plotAlignments, plotAlignments.args)
  aln_p <- aln_p + theme(legend.margin=unit(-0.5,"cm"))

  plotFreqHeatmap.args$obj = obj
  heat_p <- do.call(plotFreqHeatmap, plotFreqHeatmap.args)
  heat_p <- heat_p + theme(plot.background=element_rect(fill = "transparent",
                                                        colour = NA),
                           plot.margin = unit(c(1, 0.25, 0.5, 0), "lines"))

  arrange_args = modifyList(list(top.plot = gene_p, left.plot = aln_p,
                                 right.plot = heat_p), arrange_args)
  result <- do.call(arrangePlots, arrange_args)

  return(result)
})

#'@title Arrange plots for plotVariants:CrisprSet
#'@description Arranges 3 plots in two rows.  The vertical margins of the
#'left.plot and right.plot constrained to be equal
#'@param top.plot ggplot grob, placed on top of the figure, spanning the figure
#'width
#'@param left.plot ggplot, placed in the second row on the left
#'@param right.plot ggplot, placed in the second row on the right.
#'y-axis labels are removed.
#'@param fig.height Actual height for the figure. If not provided,
#'figure height is the sum of the row.ht.ratio (Default: NULL)
#'@param row.ht.ratio Vector specifying row height ratio (Default: c(1,6))
#'@param col.wdth.ratio  Vector specifying column width ratio (Default: c(2, 1))
#'@param left.plot.margin Unit object specifying margins of left.plot.
#'Margins of right.plot are constrained by the left.plot.
#'@return The arranged plots
arrangePlots <- function(top.plot, left.plot, right.plot, fig.height = NULL,
                      col.wdth.ratio  = c(2, 1), row.ht.ratio = c(1,6),
                      left.plot.margin = grid::unit(c(0.1,0,3,0.2), "lines")){

  # Set the size ratio of the top and bottom rows
  plot_hts <- if (is.null(fig.height)){ row.ht.ratio
  }else { fig.height/sum(row.ht.ratio)*row.ht.ratio }

  # Remove y-axis labels from right plot
  right.plot <- right.plot + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank())

  # Adjust margins of left.plot
  left.plot <- left.plot + theme(plot.margin = left.plot.margin)

  # Convert plots to grobs, lock plot heights
  p2 <- ggplot2::ggplotGrob(left.plot)
  p3 <- ggplot2::ggplotGrob(right.plot)
  p3$heights <- p2$heights

  # Return arranged plots
  return(gridExtra::grid.arrange(top.plot,
         gridExtra::arrangeGrob(p2, p3, ncol = 2, widths = col.wdth.ratio),
         nrow = 2, heights = plot_hts, newpage = FALSE))
}

#'@title Plots and annotates transcripts
#'@description Plots the gene structure, annotates this with the
#'target location
#'@param txdb A GenomicFeatures:TxDb object
#'@param target Location of target (GRanges)
#'@param target.colour Colour of box indicating targt region
#'@param target.size Thickness of box indicating target region
#'@param gene.text.size Size for figure label
#'@param panel.margin Unit object, margin size
#'@param plot.title A title for the plot.  If no plot.title is supplied,
#'the title is the list of gene ids shown (default).
#'If plot.title == FALSE, the plot will not have a title.
#'@param all.transcripts If  TRUE (default), all transcripts of genes overlapping
#'the target are shown, including transcripts that do not themselves overlap the target.
#'If FALSE, only the transcripts that overlap the target are shown.
#'@return A ggplot2 plot of the transcript structures
annotateGenePlot <- function(txdb, target, target.colour = "red",
                        target.size = 1, gene.text.size = 10,
                        panel.margin = grid::unit(c(0.1,0.1,0.1,0.1), "lines"),
                        plot.title = NULL, all.transcripts = TRUE){

  genomicfeatures <- requireNamespace("GenomicFeatures")
  stopifnot(isTRUE(genomicfeatures))

  trns <- GenomicFeatures::transcripts(txdb)
  exs <- findOverlaps(target, trns, ignore.strand = TRUE)

  if (length(exs) == 0){
    return(grid::grid.rect(gp=grid::gpar(col="white"), draw = FALSE))
  }
  # Find the genes that overlap
  genes <- AnnotationDbi::select(txdb,
              keys = as.character(mcols(trns[subjectHits(exs)])$tx_id),
              keytype = "TXID",
              columns = c("GENEID","TXNAME", "EXONSTART", "EXONEND", "TXSTRAND"))

  # Find all (possibly non-overlapping) transcripts of overlapping genes
  if (isTRUE(all.transcripts)){
    genes <- suppressWarnings(AnnotationDbi::select(txdb,
              keys = unique(genes$GENEID),
              keytype = "GENEID",
              columns = c("GENEID", "TXID", "TXNAME","EXONSTART",
              "EXONEND", "TXSTRAND")))
  }
  gene_gr <- GRanges(seqnames(target)[1], IRanges(genes$EXONSTART, genes$EXONEND),
                     strand = genes$TXSTRAND, txid = genes$TXID)

  utr5 <- GenomicFeatures::fiveUTRsByTranscript(txdb)
  utr3 <- GenomicFeatures::threeUTRsByTranscript(txdb)
  txid <- as.character(unique(genes$TXID))

  # Get UTRs matching each transcript
  all_sections <- lapply(txid, function(tid) {
    u5 <- GRangesList()
    u3 <- GRangesList()
    if (tid %in% names(utr5)){ u5 <- utr5[tid] }
    if (tid %in% names(utr3)){ u3 <- utr3[tid] }
    utrs <- unlist(c(u5, u3))
    mcols(utrs) <- NULL
    exons <- gene_gr[mcols(gene_gr)$txid == tid]
    mcols(exons) <- NULL
    split_ex <- disjoin(c(exons, utrs))
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
  tcks <- lapply(gene_spans, function(sp){
    tcks[tcks > start(sp) & tcks < end(sp)]
  })

  tck_lns <- lapply(tcks, length)
  tcks <- data.frame("tloc" = unlist(tcks),
                     "ys" = rep(1:length(tcks), lapply(tcks, length)))
  lns <- data.frame(tloc = c(start(gene_spans),end(gene_spans)),
                    ys = rep(seq_along(gene_spans),2))

  all_exs$ymax <- all_exs$ts + 0.3
  all_exs$ymin <- all_exs$ts - 0.3
  is_utr <- all_exs$type == "utr"
  all_exs$ymax[is_utr] <- all_exs$ts[is_utr] + 0.2
  all_exs$ymin[is_utr] <- all_exs$ts[is_utr] - 0.2

  target_df <- data.frame(xmin = start(target), xmax = end(target),
                          ymin = 0, ymax = ceiling(max(all_exs$ymax)))

  if (is.null(plot.title)){ plot.title <- paste(unique(genes$GENEID), sep = ";")}

  # Choose either right or left pointing arrows for the transcript plots
  strands <- unlist(lapply(all_sections, function(x) as.character(strand(x[1]))))
  strands <- rep(strands, tck_lns)
  strands[strands == "-"] <- 60
  strands[strands == "+"] <- 62
  tcks$shp <- as.integer(strands)

  p <- ggplot2::ggplot(tcks) +
    geom_point(aes_q(x = quote(tloc), y = quote(ys), group = quote(ys),
                     shape = quote(shp)), size = 2) +
    geom_line(data = lns, aes_q(x = quote(tloc), y = quote(ys),
                                group = quote(ys))) +
    scale_shape_identity()

  p <- p + geom_rect(data = all_exs, fill = "black", color = "black",
                     aes_q(xmin = quote(start), xmax = quote(end),
                           ymin = quote(ymin), ymax = quote(ymax)))

  p <- p + geom_rect(data = target_df,
                     aes_q(xmin = quote(xmin), xmax = quote(xmax),
                           ymin = quote(ymin), ymax = quote(ymax)),
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
