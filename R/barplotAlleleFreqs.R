barplotAlleleFreqs <- function(allele_counts, group = NULL, bar_colours = NULL, 
                               group_colours = NULL, legend_text_size = 10, 
                               legend_symbol_size = 1, show_percentage = TRUE, 
                               snv_label = "SNV", novar_label = "no variant"){
  
  clrs <- bar_colours
  #if (is.null(clrs)){
  #  clrs <- c("#D92120", "#E6642C", "#E68E34", "#D9AD3C", "#B5BD4C", "#7FB972",
  #            "#63AD99", "#55A1B1", "#488BC2", "#4065B1", "#413B93", "#781C81")
  #}
  
  if (is.null(clrs)){
    clrs <- c("#D92120","#E78532","#6DB388","#539EB6","#3F60AE","#781C81")
  }
  
  ac <- allele_counts
  if (!is.null(group)){
    group <- rev(group) # as ggplot plots bottom up
    if (is.null(group_colours)){
      group_colours <- c("#332288","#661100","#0072B2","#117733","#882255","#D55E00",
                         "#AA4499", "#009E73","#56B4E9","#CC79A7","#44AA99","#999933",
                         "#CC6677", "#E69F00","#88CCEE")
    }
    gp_cols <- group_colours[group]
  }
  
  #uniq <- rowSums(ac) == 1
  #ac <- ac[!uniq,, drop = FALSE]
  #ac <- rbind(ac, "Unique" = colSums(allele_counts[uniq,, drop = FALSE]))
  
  snv <- grepl(snv_label, rownames(ac))
  ac <- ac[!snv,,drop = FALSE]
  ac <- rbind(ac, "SNV" = colSums(allele_counts[snv,, drop = FALSE]))
  
  no_indel <- grepl("no variant|SNV", rownames(ac))
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
  
  var_order <- c(novar_label, snv_label, "inframe indel < 10", "inframe indel > 10",
                 "frameshift indel < 10", "frameshift indel > 10")
  
  var_labels <- c(novar_label, snv_label, "inframe indel < 10", 
                  expression("inframe indel" >= 10),
                  "frameshift indel < 10", expression("frameshift indel" >= 10))
  
  ac <- ac[intersect(var_order, rownames(ac)),]
  
  
  ### ERROR IF ONLY ONE COLUMN HERE
  af <- melt(sweep(ac, 2, colSums(ac), "/"))
  
  colnames(af) <- c("Variant", "Sample", "Percent")
  af$Variant <- factor(af$Variant, levels = var_order)                      
  
  af$Sample <- factor(af$Sample, levels = rev(unique(af$Sample)))  
  
  var_clrs <- clrs[table(af$Variant) > 0]
  
  p <- ggplot(af, aes(x = Sample, y = Percent, fill = Variant)) + 
    geom_bar(stat = "Identity", size = 10) + 
    scale_fill_manual(values = var_clrs) +  xlab(NULL) + ylab(NULL) + 
    scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
    guides(fill=guide_legend(labels = var_labels, 
                             override.aes=list(size=legend_symbol_size), nrow = 2)) + 
    theme_bw() + coord_flip() + 
    theme(legend.position = "bottom", legend.title = element_blank(), 
          legend.text = element_text(size = legend_text_size),
          plot.margin = unit(c(0.5,0.5,0.5,0),"lines"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  if (! is.null(group)){
    p <- p + theme(axis.text.y=element_text(colour= gp_cols))
    hlines <- seq_along(group)[!duplicated(group)]
    hlines <- hlines[2:length(hlines)] - 0.5
    p <- p + geom_vline(xintercept= hlines, colour = "black", size = 1)
  }
  
  #g <- ggplot(als, aes(x=1, y=1:nrow(als), label = Allele)) + geom_text(size = 4) + 
  #       geom_tile(fill = "transparent", colour = "black", size = 1) + 
  #       theme_minimal() + xlab(NULL) + ylab(NULL) + ggtitle("Variant alleles") +
  #       scale_y_continuous(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
  #       theme(axis.ticks = element_blank(), axis.text = element_blank(),
  #             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #             plot.title = element_text(vjust = -0.4, hjust = 1.1),
  #             plot.margin = unit(c(1,0,0.5,0),"lines"))
  #    
  #pgrob <- ggplotGrob(p)
  #ggrob <- ggplotGrob(g)
  #ggrob$heights <- pgrob$heights  
  #
  #grid.arrange(pgrob, ggrob, ncol = 2, widths = c(9,1))
  
  # table
  dat <- data.frame(Sample = colnames(allele_counts), 
                    Vals = c(colSums(allele_counts),colSums(allele_counts != 0)), 
                    Col = rep(c("Sequences","Alleles"), each = ncol(allele_counts)))
  
  dat$Sample <- factor(dat$Sample, level = rev(unique(af$Sample))) 
  
  q <- ggplot(dat, aes(x = Col, y = Sample, label = Vals)) + 
    geom_tile(fill = "white", colour = "black", size = 1) + geom_text() + 
    scale_y_discrete(expand = c(0,0)) + scale_x_discrete(expand = c(0,0)) +
    theme_bw() + xlab(NULL) + ylab(NULL) + 
    theme(axis.text.x = element_text(angle = 90),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = unit(c(0.25,0.25,10,0), "lines"))
  
  # q <- q + scale_fill_manual(value = "transparent")
  
  #print(q)
  
  pgrob <- ggplotGrob(p)
  ggrob <- ggplotGrob(q)
  ggrob$heights <- pgrob$heights  
  
  return(grid.arrange(pgrob, ggrob, ncol = 2, widths = c(8,2), newpage = FALSE))
  
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # g <- g + geom_text(data = counts, aes(label = Count, fill = NULL, fontface = ff), 
  #                  size = plot_text_size)
}
