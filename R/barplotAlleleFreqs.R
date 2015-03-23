#'@title Plots barplots of the spectrum of variants for a sample set
#'@author Helen Lindsay
#'@rdname barplotAlleleFreqs
#'@export
setGeneric("barplotAlleleFreqs", function(obj, ...) {
  standardGeneric("barplotAlleleFreqs")})


#'@description (signature("CrisprSet")) Groups variants by size and type 
#'and produces a barplot showing the variant spectrum for each sample.
#'Requires "VariantAnnotation" 
#'@param txdb A transcript database object
#'@rdname barplotAlleleFreqs
setMethod("barplotAlleleFreqs", signature("CrisprSet"),  
  function(obj, txdb = NULL){
    
    var_type <- obj$classifyVariantsByLoc(txdb)
    is_coding <- var_type == "coding"
    indels <- .self$cigar_freqs[is_coding,,drop = FALSE]
    
    ac <- obj$cigar_freqs
    snv.label <- obj$pars["mismatch.label"]
    novar.label <- obj$pars["match.label"]
    
    snv <- grepl(snv.label, rownames(ac))
    if (any(snv)){
      ac <- ac[!snv,,drop = FALSE]
      ac <- rbind(ac, "SNV" = colSums(obj[snv,, drop = FALSE]))
    }
    no_indel <- grepl(sprintf("%s|%s", novar.label, snv.label), rownames(ac))
    indels <- ac[!no_indel,,drop = FALSE]
    temp <- lapply(rownames(indels), function(x) strsplit(x, ",")[[1]])
    indel_grp <- rep(c(1:nrow(indels)), lapply(temp, length))
    indel_ln <- rowsum(as.numeric(gsub("^.*:([0-9]+)[DI]", "\\1", unlist(temp))), indel_grp)
    
    inframe <- indel_ln %% 3 == 0
    is_short <- indel_ln < 10  
    
    indel_grp <- rep("inframe indel < 10", nrow(indels))
    indel_grp[is_short &! inframe] <- "frameshift indel < 10"
    indel_grp[!is_short & inframe] <- "inframe indel > 10"   
    indel_grp[!is_short & !inframe] <- "frameshift indel > 10"
    
    grouped <- rowsum(indels, indel_grp)
    ac <- ac[no_indel,,drop = FALSE] 
    ac <- rbind(ac, grouped)
    
    var_order <- c(novar.label, snv.label, "inframe indel < 10", "inframe indel > 10",
                   "frameshift indel < 10", "frameshift indel > 10")
    
    var_labels <- c(novar.label, snv.label, expression("inframe indel" <= 9), 
                    "inframe indel > 10",  "frameshift indel < 9",
                    expression("frameshift indel" >= 10))
    
    names(var_labels) <- var_order
    common <- intersect(var_order, rownames(ac)) 
    ac <- ac[common,, drop=FALSE]
    af <- melt(sweep(ac, 2, colSums(ac), "/"))
    
    
})



#'@param obj
#'@param group
#'@param bar.colours
#'@param group.colours
#'@param legend.text.size
#'@param axis.text.size 
#'@param legend.symbol.size
#'@param snv.label
#'@param novar.label
#'@rdname barplotAlleleFreqs
setMethod("barplotAlleleFreqs", signature("matrix"),  
          function(obj, group = NULL, bar.colours = NULL, 
                   group.colours = NULL, legend.text.size = 10, 
                   axis.text.size= 10, legend.symbol.size = 1, 
                   snv.label = "SNV", novar.label = "no variant"){
  
  clrs <- bar.colours
  
  if (is.null(clrs)){
    clrs <- c("#D92120","#E78532","#6DB388","#539EB6","#3F60AE","#781C81")
  }
  
  ac <- obj
  if (!is.null(group)){
    group <- rev(group) # as ggplot plots bottom up
    if (is.null(group.colours)){
      group.colours <- c("#332288","#661100","#0072B2","#117733","#882255","#D55E00",
                         "#AA4499", "#009E73","#56B4E9","#CC79A7","#44AA99","#999933",
                         "#CC6677", "#E69F00","#88CCEE")
    }
    gp_cols <- group.colours[group]
  }
  
  snv <- grepl(snv.label, rownames(ac))
  if (any(snv)){
    ac <- ac[!snv,,drop = FALSE]
    ac <- rbind(ac, "SNV" = colSums(obj[snv,, drop = FALSE]))
  }
  no_indel <- grepl(sprintf("%s|%s", novar.label, snv.label), rownames(ac))
  indels <- ac[!no_indel,,drop = FALSE]
  temp <- lapply(rownames(indels), function(x) strsplit(x, ",")[[1]])
  indel_grp <- rep(c(1:nrow(indels)), lapply(temp, length))
  indel_ln <- rowsum(as.numeric(gsub("^.*:([0-9]+)[DI]", "\\1", unlist(temp))), indel_grp)
  
  inframe <- indel_ln %% 3 == 0
  is_short <- indel_ln < 10  
  
  indel_grp <- rep("inframe indel < 10", nrow(indels))
  indel_grp[is_short &! inframe] <- "frameshift indel < 10"
  indel_grp[!is_short & inframe] <- "inframe indel > 10"   
  indel_grp[!is_short & !inframe] <- "frameshift indel > 10"
  
  grouped <- rowsum(indels, indel_grp)
  ac <- ac[no_indel,,drop = FALSE] 
  ac <- rbind(ac, grouped)
  
  var_order <- c(novar.label, snv.label, "inframe indel < 10", "inframe indel > 10",
                 "frameshift indel < 10", "frameshift indel > 10")
  
  var_labels <- c(novar.label, snv.label, expression("inframe indel" <= 9), 
                  "inframe indel > 10",  "frameshift indel < 9",
                  expression("frameshift indel" >= 10))
  
  names(var_labels) <- var_order
  common <- intersect(var_order, rownames(ac)) 
  ac <- ac[common,, drop=FALSE]
  af <- melt(sweep(ac, 2, colSums(ac), "/"))
  
  colnames(af) <- c("Variant", "Sample", "Percent")
  af$Variant <- factor(af$Variant, levels = var_order)                      
  af$Sample <- factor(af$Sample, levels = rev(unique(af$Sample)))  
  var_clrs <- clrs[table(af$Variant) > 0]
  
  # barplot
  p <- ggplot(af, aes(x = Sample, y = Percent, fill = Variant)) + 
    geom_bar(stat = "Identity", size = 10) + 
    scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +   
    scale_fill_manual(values = var_clrs, labels = var_labels[common]) + 
    guides(fill=guide_legend(override.aes=list(size=legend.symbol.size), nrow = 2)) + 
    xlab(NULL) + ylab(NULL) + 
    theme_bw() + coord_flip() + 
    theme(legend.position = "bottom", legend.title = element_blank(), 
          axis.text = element_text(size = axis.text.size),
          legend.text = element_text(size = legend.text.size),
          plot.margin = unit(c(0.5,0.7,0.5,0),"lines"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  if (! is.null(group)){
    p <- p + theme(axis.text.y=element_text(colour= gp_cols))
    hlines <- seq_along(group)[!duplicated(group)]
    hlines <- hlines[2:length(hlines)] - 0.5
    p <- p + geom_vline(xintercept= hlines, colour = "black", size = 1)
  }
  # table
  dat <- data.frame(Sample = colnames(obj), 
                    Vals = c(colSums(obj),colSums(obj != 0)), 
                    Col = rep(c("Sequences","Alleles"), each = ncol(obj)))
  
  dat$Sample <- factor(dat$Sample, level = rev(unique(af$Sample))) 
  
  q <- ggplot(dat, aes(x = Col, y = Sample, label = Vals)) + 
    geom_tile(fill = "white", colour = "black", size = 1) + geom_text() + 
    scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
    theme_bw() + xlab(NULL) + ylab(NULL) + 
    theme(axis.text.x = element_text(angle = 90),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = unit(c(0.25,0.25,10,0), "lines"))
  
  pgrob <- ggplotGrob(p)
  ggrob <- ggplotGrob(q)
  ggrob$heights <- pgrob$heights  
  
  return(grid.arrange(pgrob, ggrob, ncol = 2, widths = c(8,2), newpage = FALSE))
})
