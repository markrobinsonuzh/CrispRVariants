library(GenomicAlignments)
library(Rsamtools)
library(ggplot2)
library(reshape2)
library(gridExtra)

# sangerseqR, GenomeFeatures should be in "suggests"

# Allow CrisprRun getVariants to work with a filtered variant table
# Add experiment name to CrisprSet (pars?)
# readsByPCRPrimer - could separate out searching for partial overlaps for speed
# Be consistent about target_loc / cut_site
# Give ab1ToFasta a open = "a" option to allow appending or overwriting?
# To do - add name to CrisprSet, add names to CrisprSet$cripsr_runs (easier access)
# No on target runs shouldn't stop script entirely?
# Possible bug - don't actually check the start of cigars, could have the same cigar different start?
# To do - check is.null, change to na?
# To do - design of getVarsEnsemblFormat - make fully separate from getVariants
# To do - check for consistent naming, e.g. ref versus genome
# To do - warn with plotting if multiple guides?
# Default mapping function
# TO DO - ADD THE CIGAR TO THE INSERTION TABLE IN THE CRISPR_RUN CLASS?
# To do - check that insertion_site table has "cigar" column
# Can I assume that all reads in a set have different names?
# Warning that normal cigar misleading in aln plots (or indicate start if not target_start)
# Take is.na out of plotAlignments?
# automatic layout heights and widths - absolute for gene track?
# Allow panelplot to take a gene instead of fetching all
# Remove hard coding of plots where possible
# Allow dot args in init and dispatch to appropriate function
# TO DO - make it easier to be consistent with renumbering - store option
# Multiple guide on different strands?
# pass default args to setlabels
# default param combinations - make text bigger
# store renumbered = yes / no - make sure it's consistent between functions?
# less memory if cigars weren't duplicated
# CHIMERAS: ARE OPTIONS (merge, exclude, tag) MUTUALLY EXCLUSIVE? TAG and MERGE can be together?
# Better name for classifyChimeras
# to do - why does pcdh10a fail??
# To do - check if there are more groups / insertion sites than colour (combinations)
# To do - consistency - sometimes called "target.loc", sometimes "cut.site"
# Warning if writing to non-empty file?
# Wrapper for the panel plots to ensure consistent rows?
# Check - are panelplot colours hard coded?
# getAttr function for CrisprSet lapply(cruns, function(x) x$nm, nm)?
# verbose abifToFastq


findHighCovRegions <- function(chimeras, min_cov = 1000){
  # Chimeras: GAlignments obj
  chimera_cov <- coverage(chimeras)
  sl <- slice(chimera_cov, lower = min_cov)
  slgr <- as(sl, "GRanges")
  vm <- viewMeans(sl)
  mcols(slgr)$avg <- unlist(vm)
  retunr(slgr)
}

.mergeChimeras <- function(chimeras, cigars, change_pts, unclipped, del_lens, max_overlap = 0){
  # To merge, require the indices of the chimeras in the bam, 
  # and their cigar_strings minus clipping.  Bam must have sequence available
    
  # Assumptions - not a rearrangement (the first read does not get hard clipped)
  # unclipped = clipping removed
    
  genomic_gaps <- start(chimeras[-1]) - end(chimeras[-length(chimeras)])
  
  #____________________________________
  # Change here to <= max overlap
  read_gaps <- first_aligned[-1] -  last_aligned[-length(last_aligned)] - 1
  read_gaps[!read_gaps == 0] <- sprintf("%sI", read_gaps[!read_gaps == 0])
  read_gaps[read_gaps == 0] <- ""
  #____________________________________
  
  new_g_starts <- start(chimeras[change_pts])
  new_g_ends <- end(chimeras[change_pts -1])
  
  # Keep the clipping on the end points, in case this needs to be searched for primers
  # May need to be even more specific here if end of second region is hard clipped
  new_cigars <- unclipped
  
  # First part of chimera, only clip the right:
  new_cigars[change_pts] <- gsub("(^.*M)[0-9]+[HS]","\\1", cigars[change_pts])
  
  # Last part of a chimera, only clip the left:
  new_cigars[change_pts[-1] -1] <- gsub("[0-9]+[HS](.*)", "\\1", cigars[change_pts[-1] -1])
  
  
  # HERE - DOES THIS CAUSE AN ERROR IF IT STILL IS SOFT CLIPPED?
  # If read was originally soft-clipped left but is not now, adjust the starting point
  select_start <- rep(1, length(new_cigars))  
  l_soft_clip <- grep('^[0-9]+S.*', cigars)
  select_start[l_soft_clip] <- as.numeric(gsub("(^[0-9]+)[S].*","\\1", cigars[l_soft_clip]))
  select_start[change_pts] <- 1
  select_end <- sum(width(cigarRangesAlongQuerySpace(new_cigars))) + select_start - 1


  # If one genomic location maps to multiple read locs, 
  # must trim the overlapping section to get a valid cigar string

  # OR CONSIDER IT AN INSERTION
  
  # Trim the reads
  gdels <- mcols(chimeras)$type %in% c("C:gdup", "C:rgdup")
  del_lns <- c(0, del_lns)
  rcut <- rep(0, length(del_lns))
  tocut <- del_lns <= 0
  rcut[tocut] <- -1 * del_lns[tocut] + 1
  rcut[change_pts] <- 0
  tocut <- rcut > 0 & gdels
  
  # Shift the genomic coordinates accordingly
  new_g_starts[tocut] <- new_g_starts[tocut] + rcut[tocut]
  new_g_ends[tocut] <- new_g_ends[tocut] + rcut[tocut]

  subtracted <-  as.numeric(gsub("([0-9]+)(M.*)", "\\1", new_cigars[tocut])) - rcut[tocut]
  remaining <- sprintf("%s%s", subtracted, gsub("([0-9]+)(M.*)", "\\2", new_cigars[tocut]))
  
  # When the overlap crosses multiple operations, do not merge (difficult!)
  is_positive <- !grepl("^-",subtracted)
  subtracted <- sprintf("%sI%s", rcut[tocut], remaining)
  new_cigars[tocut&is_positive] <- subtracted[is_positive] 

  sqs <- substr(mcols(chimeras)$seq, start = select_start, stop = select_end)
  sqs[change_pts - 1] <- paste0(sqs[change_pts -1], ",")


  joins <- c(sprintf("%sD%s", genomic_gaps, read_gaps), "")
  joins[change_pts[-1] -1] <- "," 
  new_cigars <- strsplit(do.call(paste0, as.list(paste0(new_cigars, joins))), ",")[[1]]

  # To do: merge the genomic duplications - find the region represented twice,
  # select only one part of read, adjust accordingly

}


getInterChimeraSeq <- function(){
  if (! "type" %in% names(cols(bam))) {
    bam <- findChimeras(bam, chimera_idxs, exclude = FALSE, tag = TRUE)
  }
  gaps <- bam[mcols(bam)$type == "C:gap"]
  # Need to get the seq minus clipping, and the gap len
    

}

classifyChimeras <- function(bam, chimera_idxs = NA, exclude = TRUE, merge = TRUE, 
                           verbose = TRUE, tag = FALSE, name = NA){    
    
    # Exclude: remove chimeras from bam file
    # Merge: join chimeras when they are long gaps
    # Tag: tag chimeras with their type
    
    # If chimera_idxs are provided, reads should be sorted by chimera name,
    # not by genomic location
    
    # Case: Aligned regions overlap
    #  1-2-3-4-5
    #      3-4-5-6-7
    #
    # Case: Inversion:
    # 1-2-3-4-5
    #            9-8-7-6 
    
    mcols(bam)$type <- NA
    
    if (exclude & tag) {
      stop("'tag' and 'exclude' are mutually exclusive.  
            Chimeras cannot be tagged if they are removed")
    }
    
    if (length(chimera_idxs) == 0) return(bam) # Length of NA = 1
    
    if (length(chimera_idxs) == 1) chimera_idxs <- findChimeras(bam) # Chimeras always >= 2       
    
    if (exclude & !merge) return(bam[-chimera_idxs])
    
    # Do all reads within a chimera map to the same chromosome?
    nms <- rle(names(bam)[chimera_idxs]) 
    nms_codes <- rep(1:length(nms$lengths), nms$lengths)
    sqs <- seqnames(bam)[chimera_idxs]
    sqs_codes <- rep(1:length(sqs@lengths), sqs@lengths)
    codes <- rle(paste(nms_codes, sqs_codes, sep = "."))
    one_chr <- rep(codes$lengths, codes$lengths) == rep(nms$lengths, nms$lengths)
    
    # And onto the same strand? (i.e. not inversion)
    strds <- strand(bam)[chimera_idxs]
    strd_rle <- rle(paste0(nms_codes, strds))
    same_strd <- rep(strd_rle$lengths, strd_rle$lengths) == rep(nms$lengths, nms$lengths)  
  
    # Are single chr chimeras gaps? (start(n+1) > end(n))
    del_lns <- start(bam)[chimera_idxs[-1]] - end(bam)[chimera_idxs[-length(chimera_idxs)]]
    is_after <- c(TRUE, del_lns > 0 )
    change_pts <- cumsum(nms$lengths) + 1 # note starts from second and includes last
    change_pts <- c(1, change_pts[1:length(change_pts) -1])
    is_after[change_pts] <- TRUE    
    codes <- rle(paste(nms_codes, is_after, sep = "."))
    has_genome_gap <- rep(codes$lengths, codes$lengths) == rep(nms$lengths, nms$lengths)
                                
    # Is the same read segment used in multiple sections of a chimera?
    # For merge-able alignments, the sum of the widths of the aligned regions of read n-1
    # should be less than or equal to the first aligned base of read n wrt the full seq
    # Note that hard-clipped regions don't appear in the cigarRanges
    cigars <- cigar(bam)[chimera_idxs]
    first_aligned <- rep(1, length(cigars))
    clipped_start <- grepl("^[0-9]+[HS]", cigars)
    
    # Note: first aligned refers to the original read, not the clipped read
    first_aligned[clipped_start] <- as.numeric(gsub("[HS].*", "", cigars[clipped_start])) + 1 
    # +1 because after the first aligned base

    # Strip the clipped bases then count the width of the remaining
    unclipped <- gsub("[0-9]+[HS]$", "", gsub("^[0-9]+[HS]", "", cigars))
    cig_ranges <- cigarRangesAlongQuerySpace(unclipped)
    last_aligned <- sum(width(cig_ranges))  
    last_aligned <- last_aligned  + first_aligned - 1
        
    # Gaps not correct: if the first part of the read maps after the 
    # second part of the read gap - however, here only care about +ve / -ve
    gaps <- c(1,  first_aligned[-1] - last_aligned[-length(last_aligned)])
    gaps[change_pts] <- 1
    gap_codes <- rle(paste(gaps > 0, nms_codes, sep = "."))
    has_read_gap <- rep(gap_codes$lengths, gap_codes$lengths) == rep(nms$lengths, nms$lengths)

    # Rearrangements: where the first aligned base of the second segment is 
    # earlier than the first 
    rearr <- first_aligned[-1] - first_aligned[-length(first_aligned)]
    rearr <- c(1, rearr)
    rearr[change_pts] <- 1
    rearr <- rearr < 0

    mergeable <- one_chr & same_strd & has_genome_gap & has_read_gap & !rearr
    
    if (verbose == TRUE){ 
      format_zero <- function(x) ifelse(is.nan(x), 0, x)
      nchm <- length(chimera_idxs)
      noc <- sum(one_chr == "TRUE")
      nss <- sum(one_chr & !same_strd == "TRUE")
      nrearr <- sum(one_chr & same_strd & rearr == "TRUE")
      ngdup <- sum(!has_genome_gap & one_chr & same_strd & !rearr == "TRUE")
      nrdup <- sum(has_genome_gap & one_chr & same_strd & !has_read_gap & !rearr == "TRUE")
      mrg <- sum(one_chr & same_strd & has_genome_gap & has_read_gap & ! rearr == "TRUE")
      if (! is.na(name)) cat(sprintf("Chimera statistics for %s:\n", name))  
      cat(sprintf(paste0("%s (%.2f%%) chimeras in %s reads\n",  
      "  %s (%.2f%%) map to the same chromosome\n", 
      "    %s (%.2f%%) map to different strands (inversions)\n",
      "    %s (%.2f%%) rearrangements (end of read maps before start)\n",
      "    %s (%.2f%%) genomic duplications (different read locs mapped to same genomic loc)\n",
      "    %s (%.2f%%) read duplications (different genomic locs mapped to same read loc)\n",
      "    %s (%.2f%%) are long gaps\n\n"),
      nchm, format_zero(nchm/length(bam)*100), length(bam),
      noc, format_zero(noc/nchm*100),
      nss, format_zero(nss/noc*100),
      nrearr, format_zero(nrearr/noc*100),
      ngdup, format_zero(ngdup/noc*100),
      nrdup, format_zero(nrdup/noc*100),
      mrg, format_zero(mrg/noc*100))) 
    }
      
   if (tag == TRUE){
     mcols(bam)$type <- "NA"
     mcols(bam)$type[chimera_idxs] <- "C"
     mcols(bam)$type[chimera_idxs[! one_chr]] <- "C:multichr"
     mcols(bam)$type[chimera_idxs[one_chr & ! same_strd]] <- "C:inv"
     mcols(bam)$type[chimera_idxs[one_chr & same_strd & ! has_genome_gap]] <- "C:gdup"
     mcols(bam)$type[chimera_idxs[one_chr & same_strd & ! has_read_gap]] <- "C:rdup"
     rg_dup <- chimera_idxs[one_chr & same_strd & ! has_genome_gap & ! has_read_gap]
     mcols(bam)$type[rg_dup] <- "C:rgdup"
     mcols(bam)$type[rearr & one_chr & same_strd] <- "C:rearr"
     mcols(bam)$type[chimera_idxs[mergeable]] <- "C:gap"
     return(bam)
   }
    
    ######
    # To do - merge these
    
    bam <- bam[-chimera_idxs]    
    return(bam)
}



plotAlleleFreqs <- function(allele_freqs, size = pt_size){
  alleles <- colSums(allele_freqs != 0)  
  als <- data.frame(Allele = alleles, Sample =names(alleles))
  p <- ggplot(als, aes(x=Sample, y = Allele, group = Group, colour = Group)) 
      + geom_line() + geom_point(size = pt_size) + theme_bw()

}

cigarFrequencyHeatmap <- function(cigar_freqs, as_percent = TRUE, x_size = 16, y_size = 16, 
                           x_axis_title = NULL, x_angle = 90, annotate_counts = TRUE,
                           col_sums = TRUE, colours = "yellowred", legend_text_size = 16,
                           plot_text_size = 8, group = NULL, group_colours = NULL){

  # legend_text_size applies to title and text, alternatively can be changed 
  # in finished plot using "theme"
      
  cig_freqs <- cigar_freqs
  if (col_sums == TRUE){
    # Make space for totals to be added
    cig_freqs <- rbind(Total = rep(NA, ncol(cig_freqs)), cig_freqs)
  }
  
  # If a sample group is supplied, reorder the columns of counts  
  if (!is.null(group)){
    group <- as.factor(group)  
    cig_freqs <- cig_freqs[,order(group), drop = FALSE]

    if (is.null(group_colours)){
              
      clrs <- c("#332288","#661100","#117733","#882255","#D55E00", 
                "#0072B2","#AA4499","#009E73","#56B4E9","#CC79A7",
                "#44AA99","#999933","#CC6677", "#E69F00","#88CCEE")
      clrs <- clrs[group]
    } else { 
      clrs <- group_colours[group] 
    }   
  }
  
  counts <- melt(cig_freqs)  
  colnames(counts) <- c("Cigar", "Sample","Count")
  counts$Cigar <- factor(counts$Cigar, levels = rev(levels(counts$Cigar)))


  if (as_percent == TRUE){
    m <- apply(cigar_freqs, 2, function(x) x/sum(x))  
    if (col_sums == TRUE){
      m <- rbind(Total = rep(NA, ncol(m)), m)
    }
    m <- melt(m)
    colnames(m) <- c("Cigar", "Sample","Percentage")  
    m$Cigar <- factor(m$Cigar, levels = rev(levels(m$Cigar)))
    g <- ggplot(m, aes(x = Sample, y = Cigar, fill = Percentage)) + geom_tile()
  }
  else{    
    g <- ggplot(counts, aes(x = Sample, y = Cigar, fill = Count)) + geom_tile()
  }
  
  if (annotate_counts == TRUE){
    counts$ff <- "plain"
    if (col_sums == TRUE){ 
      idxs <- which(is.na(counts$Count))
      counts$Count[idxs] <- colSums(cigar_freqs)
      counts$ff[idxs] <- "bold"
      
      xranges <- ggplot_build(g)$panel$ranges[[1]]$x.range  
      yranges <- ggplot_build(g)$panel$ranges[[1]]$y.range  
      box_coords <- data.frame(xmin = min(xranges), xmax = max(xranges), 
                               ymin = nrow(cig_freqs) -0.5, ymax = max(yranges))
                               
      g <- g + geom_rect(data=box_coords, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax,
                        ymax = ymax, x = NULL, y = NULL),
                         color = "black", size = 1, fill = "transparent")       
    }
    g <- g + geom_text(data = counts, aes(label = Count, fill = NULL, fontface = ff), 
                       size = plot_text_size)
  }
 
  if (colours == "yellowred"){     
    hmcols<-colorRampPalette(c("white","gold","orange","orangered","red", "darkred"))(50) 
    g <- g + scale_fill_gradientn(colours = hmcols, na.value = "white") 
  }   
  
  g <- g + ylab(NULL) + xlab(x_axis_title) + theme_bw() +
       theme(axis.text.x = element_text(size = x_size, angle = x_angle, hjust = 1),
             axis.text.y = element_text(size = y_size),
             legend.text = element_text(size = legend_text_size),
             legend.title = element_text(size = legend_text_size),
             legend.key.height = unit(5, "lines"))         
             
  if (! is.null(group)){
    g <- g + theme(axis.text.x=element_text(colour= clrs))
  }
  
  return(g)
}


getCodonFrame <- function(txdb, target_chr, target_start, target_end){
  # Target_chr must match txdb
  require(VariantAnnotation)
  refLocsToLocalLocs(GRanges(target_chr, IRanges(target_start, target_end)), txdb)
  
  # check that all transcripts have the same frame
  
}

nucleotideToAA <- function(seqs, txdb, target_start, target_end){
    
}

panelPlot <- function(txdb, target_chr, target_start, target_end, aln_p, heat_p, 
                      col_widths = c(2, 1), fig_height = NULL, row_ht_ratio = c(1,6), 
                      panel_margin = unit(c(0.25,0.25,0.25,0.25), "lines"), gene_text_size = 20, 
                      aln_p_margin = unit(c(0.25,0,10,0.5), "lines")){

  # Plot_margin must be a unit object, and will be applied to the alignment plot.  
  # The vertical margins of the heatmap are constrained to equal those of the alignment plot

  aln_hts <- if (is.null(fig_height)){ row_ht_ratio 
             }else { fig_height/sum(row_ht_ratio)*row_ht_ratio } 

  # lock the plot area heights of the alignment and heatmap
  heat_p <- heat_p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  
  #aln_p <- aln_p + annotation_custom(grob = linesGrob(), xmin = -Inf, xmax = 5, 
  #     ymin = -Inf, ymax = 10)             
             
  aln_p <- aln_p + theme(plot.margin = aln_p_margin)            
  p2 <- ggplotGrob(aln_p)
  #p2$layout$clip[p2$layout$name=="panel"] <- "off"
  
  p3 <- ggplotGrob(heat_p)
  #p3$layout$clip[p3$layout$name=="panel"] <- "off"
  #p2$heights <- p3$heights  
  
  p3$heights <- p2$heights  
  
  
  # Make the gene plot
  require(ggbio)
  genes <- genes(txdb)
  target <- GRanges(target_chr, IRanges(target_start, target_end))
  wh <- genes[findOverlaps(genes, target)@queryHits]  
    
  # ACCOUNT FOR GUIDES NOT IN A GENE?
  if (length(wh) == 0){
    p1 <- grob()
  } else{
    gene_id <- do.call(paste, mcols(wh)$gene_id)
    p1 <- autoplot(txdb, wh, label = FALSE)
    #Pull off the y limits from the transcript plot
    yranges <- ggplot_build(p1)$panel$ranges[[1]]$y.range 
    xranges <- ggplot_build(p1)$panel$ranges[[1]]$x.range
    target_df <- data.frame(xmin = start(target), xmax = end(target), 
                           ymin = yranges[1], ymax = yranges[2])
                        
    p1 <-  p1 + geom_rect(data = target_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                 colour = "red", fill = NA, size = 1) + theme_minimal() + ggtitle(gene_id) + 
                 theme(axis.text.x = element_text(size = gene_text_size), 
                 text = element_text(size = gene_text_size),
                 panel.margin = panel_margin) 
     # Add some space at the bottom of the panel 
    #p1 <- p1 + theme(plot.margin = unit(c(1,1,10,1), "lines"))
    p1 <- as.list(attributes(p1))$ggplot
    #p1$layout$clip[p1$layout$name=="panel"] <- "off"
    p1 <- ggplotGrob(p1)
   
  }
   
  # newpage heatmap means names clip again
  
  # nested arrange grob?
  return(grid.arrange(p1, arrangeGrob(p2, p3, ncol = 2, widths = col_widths), 
         nrow = 2, heights = aln_hts, newpage = FALSE))
  
  #return(arrangeGrob(p1, arrangeGrob(p2, p3, ncol = 2, widths = col_widths), 
  #       nrow = 2, heights = c(gen_ht, aln_ht), newpage = FALSE))
  
}


