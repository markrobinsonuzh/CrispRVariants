#'@title Plots barplots of the spectrum of variants for a sample set
#'@author Helen Lindsay
#'@description For signature "matrix", this function optionally does a very
#'naive classification of variants by size.  Frameshift variant combinations are
#'those whose sum is not divisible by three.  Intron boundaries are *NOT* considered,
#'use with caution!
#'For signature "CrisprSet", the function uses the VariantAnnotation package to
#'localize variant alleles with respect to annotated transcripts.  Variants are 
#'annotated as "coding" when they are coding in any transcript.
#'@rdname barplotAlleleFreqs
#'@param obj The object to be plotted
#'@param ... additional arguments
#'@return A ggplot2 barplot of the allele distribution and optionally a table
#'of allele counts
#'@export
setGeneric("barplotAlleleFreqs", function(obj, ...) {
  standardGeneric("barplotAlleleFreqs")})


#'@description (signature("CrisprSet")) Groups variants by size and type
#'and produces a barplot showing the variant spectrum for each sample.
#'Accepts all arguments accepted by barplotAlleleFreqs for signature("matrix").
#'Requires package "VariantAnnotation"
#'@param txdb A transcript database object
#'@param min.freq Include variants with at frequency least min.freq in at least
#'one sample.  (Default: 0, i.e. no cutoff)
#'@param include.chimeras Should chimeric reads be included in results?
#'(Default: TRUE)
#'@param group A grouping variable for colouring the sample names.  
#'(Default: NULL)
#'@rdname barplotAlleleFreqs
setMethod("barplotAlleleFreqs", signature("CrisprSet"),
  function(obj, ..., txdb, min.freq = 0, include.chimeras = TRUE, group = NULL){

    # Potential improvements:
    # Do filtering before variant location
    
    clrs <- c("#D92120","#E6642C","#488BC2","#4065B1","#413B93",
              "#781C81","#E68E34","#D9AD3C","#B5BD4C","#7FB972",
              "#63AD99","#55A1B1","Gray")
    
    var_order <- c(obj$pars["match_label"], obj$pars["mismatch_label"],
                   "inframe indel < 10", "inframe indel > 10",
                   "frameshift indel < 10", "frameshift indel > 10",
                   "spliceSite","intron","fiveUTR",
                   "threeUTR","promoter", "intergenic", "Chimeric")
    
    var_labels <- c(obj$pars["match_label"], obj$pars["mismatch_label"],
                    expression("inframe indel" <= 9),
                    "inframe indel > 10",  "frameshift indel < 9",
                    expression("frameshift indel" >= 10),
                    "splice site","intronic","5' UTR",
                    "3' UTR","promoter", "intergenic", "Chimeric")
    
    var_type <- obj$classifyVariantsByLoc(txdb)
    classification <- obj$classifyCodingBySize(var_type)
    classification["Other"] <- "Chimeric"
    
    # Remove variants below the min.freq
    ac <- obj$.getFilteredCigarTable(min.freq = min.freq,
                                     include.chimeras = include.chimeras)
    rownames(ac) <- classification[rownames(ac)]
    vls <- factor(rownames(ac), levels = var_order)
    ac <- rowsum(ac, vls)
    clrs <- clrs[var_order %in% rownames(ac)]
    barplotAlleleFreqs(ac, bar.colours = clrs, group = group,
                       classify = FALSE, ...)
})


#'@description signature("matrix") Accepts a matrix of allele counts,
#'with rownames being alleles and column names samples.
#'@param group A grouping factor for the columns in obj.  Columns in the
#'same group will be displayed in the same text colour (Default: NULL)
#'@param bar.colours Colours for the categories in the barplot.
#'Colours must be provided if there are more than 6 different categories.
#'@param group.colours Colours for the text labels for the experimental groups
#'A set of 15 different colours is provided.
#'@param legend.text.size The size of the legend text, in points.
#'@param axis.text.size The size of the axis text, in points
#'@param legend.symbol.size  The size of the symbols in the legend
#'@param snv.label The row label for single nucleotide variants
#'@param novar.label The row label for non-variant sequences
#'@param include.table Should a table of allele (variant combination)
#'counts and total sequences be plotted? (Default: TRUE)
#'@param chimera.label The row label for chimeric (non-linearly aligned)
#'variant alleles
#'@param classify If TRUE, performs a naive classification by size 
#'(Default:TRUE)
#'@rdname barplotAlleleFreqs
#'@examples
#'data("gol_clutch1")
#'barplotAlleleFreqs(variantCounts(gol))
#'
#'# Just show the barplot without the counts table:
#'barplotAlleleFreqs(variantCounts(gol), include.table = FALSE)
setMethod("barplotAlleleFreqs", signature("matrix"),
          function(obj, group = NULL, bar.colours = NULL,
                   group.colours = NULL, legend.text.size = 10,
                   axis.text.size = 10, legend.symbol.size = 1,
                   snv.label = "SNV", novar.label = "no variant",
                   chimera.label = "Other", include.table = TRUE,
                   classify = TRUE){

  clrs <- bar.colours

  if (is.null(clrs)){
    # Colour-blind safe rainbow palette
    clrs <- c("#D92120","#E78532","#DFA53A", "#6DB388","#539EB6",
              "#3F60AE","#781C81")
  }

  ac <- obj
  if (!is.null(group)){
    group <- rev(group) # as ggplot plots bottom up
    if (is.null(group.colours)){
      # A range of dark colours chosen for readability and distinctiveness
      group.colours <- c("#332288","#661100","#0072B2","#117733","#882255",
                         "#D55E00","#AA4499","#009E73","#56B4E9","#CC79A7",
                         "#44AA99","#999933","#CC6677", "#E69F00","#88CCEE")
      }
    if (length(levels(group)) > length(group.colours)){
      stop("Too many groups for default colours.  Supply a vector of
           group.colours with length equal to the number of groups")
    }
    gp_cols <- group.colours[group]
  }
  
  if (isTRUE(classify)){
    temp <- classifyBySize(ac, snv.label, novar.label, chimera.label, clrs)
    ac <- temp$ac
    var_labels <- temp$var_labels
    var_clrs <- temp$var_clrs
  } else {
    var_labels <- rownames(ac)
    var_clrs <- clrs
  }

  af <- reshape2::melt(sweep(ac, 2, colSums(ac), "/"))

  colnames(af) <- c("Variant", "Sample", "Percent")
  af$Variant <- factor(af$Variant, levels = rownames(ac))
  af$Sample <- factor(af$Sample, levels = rev(unique(af$Sample)))

  # barplot
  p <- ggplot(af, aes_q(x = quote(Sample), y = quote(Percent),
                        fill = quote(Variant))) +
    geom_bar(stat = "Identity", size = 10) +
    scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
    scale_fill_manual(values = var_clrs, labels = var_labels) +
    guides(fill=guide_legend(override.aes=list(size=legend.symbol.size),
                             ncol = 3, label.hjust = 0)) +
    xlab(NULL) + ylab(NULL) +
    theme_bw() + coord_flip() +
    theme(legend.position = "bottom", legend.title = element_blank(),
          axis.text = element_text(size = axis.text.size),
          legend.text = element_text(size = legend.text.size),
          plot.margin = grid::unit(c(0.5,0.7,0.5,0),"lines"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  if (! is.null(group)){
    p <- p + theme(axis.text.y=element_text(colour= gp_cols))
    hlines <- seq_along(group)[!duplicated(group)]
    hlines <- hlines[2:length(hlines)] - 0.5
    p <- p + geom_vline(xintercept= hlines, colour = "black", size = 1)
  }
  # table
  if (include.table == FALSE){ return(p) }

  dat <- data.frame(Sample = colnames(obj),
                    Vals = c(colSums(obj),colSums(obj != 0)),
                    Col = rep(c("Sequences","Alleles"), each = ncol(obj)))

  dat$Sample <- factor(dat$Sample, levels = rev(unique(af$Sample)))

  q <- ggplot(dat, aes_q(x = quote(Col), y = quote(Sample),
                         label = quote(Vals))) +
    geom_tile(fill = "white", colour = "black", size = 1) + geom_text() +
    scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
    theme_bw() + xlab(NULL) + ylab(NULL) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = grid::unit(c(0.25,0.25,10,0), "lines"))

  pgrob <- ggplot2::ggplotGrob(p)
  ggrob <- ggplot2::ggplotGrob(q)
  ggrob$heights <- pgrob$heights

  return(gridExtra::grid.arrange(pgrob, ggrob, ncol = 2,
                                 widths = c(8,2), newpage = FALSE))
})


aggregateSNVs <- function(ac, snv.label){
  snv <- grepl(snv.label, rownames(ac))
  if (any(snv)){
    temp <- ac[!snv,,drop = FALSE]
    ac <- rbind(temp, "SNV" = colSums(ac[snv,, drop = FALSE]))
  }
  ac
}


classifyBySize <- function(ac, snv.label, novar.label, chimera.label, clrs){
  ac <- aggregateSNVs(ac, snv.label)
  no_indel <- grepl(sprintf("%s|%s|%s", novar.label, snv.label, chimera.label),
                    rownames(ac))
  
  # Classify indel variants by size
  if (any(!no_indel)){
    indels <- ac[!no_indel,,drop = FALSE]
    temp <- strsplit(rownames(indels), ",")
    indel_grp <- rep(c(1:nrow(indels)), lapply(temp, length))
    indel_ln <- rowsum(as.numeric(gsub("^.*:([0-9]+)[DI]", "\\1", unlist(temp))),
                       indel_grp)
    
    inframe <- indel_ln %% 3 == 0
    is_short <- indel_ln < 10
    
    indel_grp <- rep("inframe indel < 10", nrow(indels))
    indel_grp[is_short &! inframe] <- "frameshift indel < 10"
    indel_grp[!is_short & inframe] <- "inframe indel > 10"
    indel_grp[!is_short & !inframe] <- "frameshift indel > 10"
    
    grouped <- rowsum(indels, indel_grp)
    ac <- ac[no_indel,,drop = FALSE]
    ac <- rbind(ac, grouped)
  }
  
  var_order <- c(novar.label, snv.label, chimera.label,
                 "inframe indel < 10", "inframe indel > 10",
                 "frameshift indel < 10", "frameshift indel > 10")
  
  var_labels <- c(novar.label, snv.label, chimera.label,
                  expression("inframe indel" <= 9),
                  "inframe indel > 10",  "frameshift indel < 9",
                  expression("frameshift indel" >= 10))
  
  names(var_labels) <- var_order
  common <- intersect(var_order, rownames(ac))
  var_labels <- var_labels[common]
  ac <- ac[common,, drop=FALSE]
  var_clrs <- clrs[factor(common, levels = var_order)]
  list(ac = ac, var_labels = var_labels, var_clrs = var_clrs)
}
