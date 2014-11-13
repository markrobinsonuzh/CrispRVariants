library(GenomicAlignments)
library(Rsamtools)
library(ggplot2)
library(reshape2)
library(gridExtra)

# sangerseqR, GenomeFeatures should be in "sugests"

# Give ab1ToFasta a open = "a" option to allow appending or overwriting?
# CrisprSet needs to store the guide location, cut site
# To do - add name to CrisprSet, add names to CrisprSet$cripsr_runs (easier access)
# No on target runs shouldn't stop script entirely?
# Possible bug - don't actually check the start of cigars, could have the same cigar different start?
# To do - store target in CrisprSet (lose flexibility to specify a smaller set)
# To do - check is.null, change to na?
# To do - design of getVarsEnsemblFormat - make fully separate from getVariants
# To do - check for consistent naming, e.g. ref versus genome
# To do - deal with paired 
# To do - warn with plotting if multiple guides?
# Default mapping function
# TO DO - ADD THE CIGAR TO THE INSERTION TABLE IN THE CRISPR_RUN CLASS?
# To do - check that insertion_site table has "cigar" column
# To do - indel realignment around target loc
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
# To do - consistency - sometimes called "target_loc", sometimes "cut_site"
# Warning if writing to non-empty file?
# Wrapper for the panel plots to ensure consistent rows?
# Check - are panelplot colours hard coded?

# Show method for CrisprRun, CrisprMultiplex
# Filter by mapq for CrisprRun and argument for CrisprSet, as in CrisprMultiplex
# getAttr function for CrisprSet lapply(cruns, function(x) x$nm, nm)?
# verbose abifToFastq

amino_colours <- matrix(c("H", "#5555ff", "#8282D2", "#7070FF",
                          "K", "#ffcc77", "#145AFF", "#4747B8",
                          "R", "#ffcc77", "#145AFF", "#00007C",
                          "D", "#55bb33", "#E60A0A", "#A00042",
                          "E", "#55bb33", "#E60A0A", "#660000",
                          "S", "#ff4455", "#FA9600", "#FF4C4C",
                          "T", "#ff4455", "#FA9600", "#A00042",
                          "N", "#55bb33", "#00DCDC", "#FF7C70",
                          "Q", "#55bb33", "#00DCDC", "#FF4C4C",
                          "A", "#77dd88", "#C8C8C8", "#8CFF8C",
                          "V", "#66bbff", "#0F820F", "#FFFFFF",
                          "L", "#66bbff", "#0F820F", "#455E45",
                          "I", "#66bbff", "#0F820F", "#004C00",
                          "M", "#66bbff", "#E6E600", "#B8A042",
                          "F", "#9999ff", "#3232AA", "#534C42",
                          "Y", "#9999ff", "#3232AA", "#B8A042",
                          "W", "#9999ff", "#B45AB4", "#534C42",
                          "P", "#eeaaaa", "#DC9682", "#534C42",
                          "G", "#77dd88", "#EBEBEB", "#FFFFFF",
                          "C", "#99ee66", "#E6E600", "#FFFF70"), 
                          byrow = TRUE, ncol = 4)
                         
colnames(amino_colours) <- c("AA", "MAEditor","Amino","Shapely")
row.names(amino_colours) <- amino_colours[,"AA"]
#"B",
#"Z",
#"X")


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

plotAlignments <- function(ref, alns, ins_sites, pam_loc = NA, show_plot = FALSE, 
                           target_loc = 18, pam_start = NA, pam_end = NA, 
                           ins_size = 6, legend_cols = 3, xlab = NULL, xtick_labs = NULL,
                           xtick_breaks = NULL, plot_text_size = 8, axis_text_size = 16, 
                           legend_text_size = 16, highlight_guide=TRUE, guide_loc = NULL){
  
  # ref: the reference sequence
  # alns: named vector of aligned sequences, with insertions removed
  # ins_sites:  table of insertion_sites, must include cols named "start", "cigar" and "seq"
  # pam_loc: location of PAM with respect to the target site
  # Insertion locations are determined by matching ins_sites$cigar with names(alns)
  # All characters other than ACTG are labelled N
    
  # Reverse alignment order, as ggplot geom_tile plots bottom up  
  aln_chrs <- strsplit(c(rev(alns), Reference = as.character(ref)), "")
  
  # Test that all alignments and reference have the same length
  if (! length(unique(lapply(aln_chrs, length))) == 1){
    stop("The reference sequence and the alignments must all have the same length")
  }
    
  temp <- t(as.data.frame(aln_chrs))
  rownames(temp) <- names(aln_chrs)
  m <- melt(temp)
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
  
  # Colours and shapes for the insertion markers and tiles
  shps <- c(21,23,25) 
  clrs <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
            "#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100",
            "#CC6677","#882255", "#AA4499")


  # Plot aligned sequences
  tile_height <- 0.55
     
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
  
  # Add line for the cut site
  p <- p + geom_vline(xintercept= target_loc + 0.5, colour = "black", size = 2)# linetype = "dashed",
                                     
  # Make a data frame of insertion locations 
  if (nrow(ins_sites) > 0){
    
    ins_ord <- match(ins_sites$cigar, names(aln_chrs))    
  
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
  ymin = length(names(aln_chrs)) - (tile_height / 2 + 0.25)
  ymax = length(names(aln_chrs)) + (tile_height / 2 + 0.25)
  
  if (is.null(guide_loc)){
    xmin <- target_loc - 17.5
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
  if (! is.na(pam_start)){    
    if (is.na(pam_end)){
      pam_end <- pam_start + 2
    }
    pam_df <- data.frame(xmin = pam_start - 0.5, xmax = pam_end - 0.5,
                         ymin = length(names(aln_chrs)) - (tile_height / 2 + 0.2),
                         ymax = length(names(aln_chrs)) + (tile_height / 2 + 0.2))
    p <- p + geom_rect(data=pam_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax,
                        ymax = ymax, x = NULL, y = NULL),
                         color = "black", size = 3, fill = "transparent") 
    p <- p + annotation_custom(grob = textGrob("PAM", gp = gpar(cex = 3)), 
               xmin = pam_df$xmin, xmax = pam_df$xmax, 
               ymin = pam_df$ymin + 1, ymax = pam_df$ymax + 1)
  }
      

  if (show_plot == TRUE){
    print(p)
  }
  return(p)
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


#_______________________________________________________________________________________




#_______________________________________________________________________________________

CrisprRun = setRefClass(
  Class = "CrisprRun",
  fields = c("alns",
             "name",
             "query_ranges",
             "ref_ranges",
             "genome_ranges",
             "cigar_ops",
             "off_target",
             "long_dels",
             "insertions",
             "ins_key", 
             ".read_types",
             "cigar_labels")
)

CrisprRun$methods(
  initialize = function(bam, target_chr, target_start, target_end, rc = FALSE, 
                        name = NULL, count_unmapped = TRUE, merge_chimeras = FALSE, 
                        tag_chimeras = TRUE, exclude_chimeras = FALSE, 
                        exclude_ranges = GRanges(), exclude_names = NA, verbose = TRUE){
    
    # bam can either be GAlignments or the filename to read
    # Cigars don't work as is for plotting, so set cigar_label attribute                    
    
    .read_types <<- list()
    
    name <<- ifelse(is.null(name), bam, name)

    if (verbose == TRUE) cat(sprintf("Initialising %s\n", .self$name))

    if (class(bam) == "character"){
    
      galns <- .readBam(bam, target_chr, target_start, target_end, merge_chimeras, 
                  tag_chimeras, exclude_chimeras, exclude_ranges, exclude_names, verbose)
    
      if (count_unmapped == TRUE){
        param <- ScanBamParam(what = c("qname"), flag=scanBamFlag(isUnmappedQuery=TRUE))
        .read_types[scanBam(bam, param = param)[[1]]$qname] <<- "unmapped"
      }
        
      off_target <<- galns[["off_target"]]
      on_target <- galns[["on_target"]]
        
      if (length(on_target) == 0) return()
    
    } else if (class(bam) == "GAlignments"){
      on_target <- bam
    
    } else {
      stop("bam can only be a filename or a GAlignments object")
    }
    ref_ranges <<- cigarRangesAlongReferenceSpace(cigar(on_target))
    query_ranges <<- cigarRangesAlongQuerySpace(cigar(on_target))
    
    cigar_ops <<- CharacterList(explodeCigarOps(cigar(on_target)))
    .self$readsToTarget(on_target, target_start, target_end, rc, verbose)
    .self$getInsertionSeqs()
  },
  
  show = function(){
    print(c(class(.self), sprintf("CrisprRun object named %s, with %s on target alignments.", 
                .self$name, length(.self$alns)), .self$alns))
  },
  
  .readBam = function(bam_fname, target_chr, target_start, target_end, merge_chimeras,
                     tag_chimeras, exclude_chimeras, exclude_ranges = GRanges(), 
                     exclude_names = NA, verbose = TRUE){
                     
    # This function reads a bam file and classifies reads as on_target if they 
    # completely span the target, region, and off_target otherwise. 
    #
    # Note that unmapped reads are discarded by readGAlignments, 
    # and are not counted as "off-target" reads
    #
    # "exclude_ranges" should be a GRanges object of regions that should not be counted,
    #  e.g. primer or cloning vector sequences that have a match in the genome   
    #
    # If merge_chimeras = TRUE, chimeric reads (alignments with a large gap, split 
    # over two lines) are merged if both sequences are from the same chromosome, 
    # do not overlap, and are aligned to the same strand.
    # It is assumed that sequences with two alignments are chimeras, not alternate mappings
    # 
    # All alignments are read in, and off-target alignments are annotated    
    vb <- verbose
  
    param <- ScanBamParam(what = c("seq"))
    bam <- readGAlignments(bam_fname, param = param, use.names = TRUE)
    
    #.read_types <<- sapply(unique(names(bam)), function(x) NA) 
    
    #Exclude reads
    temp <- excludeFromBam(bam, exclude_ranges, exclude_names)    
    .read_types[setdiff(names(bam), names(temp))] <<- "excluded"
    if (vb == TRUE){
       original <- length(bam)
       cat(sprintf("Read %s alignments, excluded %s\n", original, original - length(temp)))
    }
    bam <- temp         
    
    # Tag, merge, and/or exclude chimeric sequences
    if (length(bam) == 0) return(list("on_target" = NULL, "off_target" = NULL))
   
        
    if (any(c(merge_chimeras, exclude_chimeras, tag_chimeras)) == TRUE){      
      if (tag_chimeras == TRUE){
        chimera_idxs <- findChimeras(bam)
        bam <- classifyChimeras(bam, chimera_idxs, exclude = FALSE, merge = FALSE, 
                              tag = tag_chimeras, verbose = vb, name = name)
        
        if (merge_chimeras == TRUE){
        
        }
        if (exclude_chimeras == TRUE){
          chimera_nms <- names(bam[chimera_idxs])
          chimera_types <- mcols(bam[chimera_idxs])$type
          
          # ADD UNIQUE COMBINATIONS OF THESE TO READ_TYPES ARG
        }
        
      } else {
        bam <- classifyChimeras(bam, exclude = exclude_chimeras, merge = merge_chimeras, 
                              tag = tag_chimeras, verbose = vb, name = name)
      }
      if (verbose == TRUE){
        cat(sprintf("%s reads after filtering chimeras\n", length(bam)))
      }
      if (length(bam) == 0) return(list("on_target" = NULL, "off_target" = NULL))
    }
    
    # Find reads that partially cover the target site 
    partially_on_target <- factor(as.character((start(bam) <= target_start | end(bam) >= target_end) & 
                     seqnames(bam) == target_chr), levels = c("TRUE", "FALSE"))
    
    partial_names <- split(names(bam), partially_on_target)    
    names(partial_names) <- c("partial", "off_target")  
       
    # Split reads into fully on-target and off-target    

    is_on_target <- factor(as.character(start(bam) <= target_start & end(bam) >= target_end & 
                     seqnames(bam) == target_chr), levels = c("TRUE", "FALSE"))
    result <- split(bam, is_on_target)
    names(result) <- c("on_target", "off_target")
    
    if (vb == TRUE){
        cat(sprintf("  %s on-target, %s off-target or partially on-target\n",
                       length(result$on_target), length(result$off_target))) }
                        
    # Note that a read can currently be classified as "excluded" if it's a chimera
    # and half is excluded
    
    remaining <- names(which(is.na(.read_types) | .read_types == "excluded"))
    .read_types[remaining[remaining %in% names(result$off_target)]] <<- "off_target"
    .read_types[remaining[remaining %in% partial_names$partial]] <<- "partially_on_target"
    .read_types[remaining[remaining %in% names(result$on_target)]] <<- "on_target"
    result
  },   
  
  readsToTarget = function(alns, target_start, target_end, rc, verbose = TRUE){
    # Narrow the reads and cigar strings to the target region, set these attributes
    
    # Narrowing example:
    # 3-4-5-6-7-8-9-10 Read
    #     5-6-7-8      Target sequence
    # target_start 5 - (read_start 3 - 1) = index 3
    # target_end 8 - target_start 5 + cigstart 3 = index 6
    # Note that alignments and cigars are reversed if strand is -ve,
    #  but start is still genomic
    
    if (verbose == TRUE) cat("finding deletions \n")
    locs <- .self$findDeletions(target_start, target_end, alns)
    
    if (verbose == TRUE) cat("narrowing\n")
    # Record how many bases before first aligned read    
    clip_starts <- rep(0, length(alns))
    is_clipped <- grepl("^[0-9]+[HS]", cigar(alns))
    clip_starts[is_clipped] <- unlist(lapply(width(.self$query_ranges[is_clipped]), '[[', 1))
            
    temp <- cigarNarrow(cigar(alns), locs$starts, locs$ends)  
    new_starts <-  attr(temp, "rshift") + 1 + clip_starts 
    # + 1 as rshift is number removed not starting point
    # new_starts are read offsets wrt current genomic starts.
    
    # Account for insertions and deletions prior to the target site
    # Need to adjust for deletions as new_starts derived wrt reference
    prior <- start(.self$ref_ranges) <= locs$start
    prior_ins <- prior & .self$cigar_ops == "I"
    prior_del <- prior & .self$cigar_ops == "D" 
    new_starts <- new_starts + sum(width(.self$query_ranges)[prior_ins])
    new_starts <- new_starts - sum(width(.self$ref_ranges)[prior_del])
    
    # Set cigars and update ranges attributes
    cigs <- as.character(temp)
    query_ranges <<- cigarRangesAlongQuerySpace(cigs) 
    shift_starts <- start(alns) + attr(temp, "rshift") -1
    genome_ranges <<- shift(cigarRangesAlongReferenceSpace(cigs), shift_starts)
        
    seq_lens <- sum(width(.self$query_ranges))  
    seqs <- subseq(mcols(alns)$seq, start = new_starts, width = seq_lens)
    genome_start <- start(alns) + attr(temp, "rshift")
    
    if (rc == TRUE){
       if (verbose == TRUE) print("reversing narrowed alignments")
       cigs <- unname(sapply(cigs, reverseCigar))
       query_ranges <<- cigarRangesAlongQuerySpace(cigs)
       # Note that even if strand is -ve, it is displayed wrt reference strand in bam
       seqs <- reverseComplement(seqs)
       genome_ranges <<- IRangesList(lapply(.self$genome_ranges, rev))
       
       # Add difference between target start and actual start, subtract difference at end
       lshift <- target_start - (start(alns) - 1) - locs$starts
       rshift <- target_end - target_start - ( locs$ends - locs$starts - lshift)           
       genome_start <- as.integer(genome_start + lshift + rshift)
    }
   
    ref_ranges <<- cigarRangesAlongReferenceSpace(cigs)
    cigar_ops <<- CharacterList(explodeCigarOps(cigs))
    alns <<- GAlignments(seqnames = seqnames(alns), pos= genome_start, cigar = cigs, 
                names = names(alns), strand = strand(alns), seqlengths = seqlengths(alns),
                seq = seqs, type = mcols(alns)$type) 
   },
  
  countIndels = function(){
    return(sum(any(.self$cigar_ops %in% c("I", "D"))))
  },
  
  indelPercent = function(){
    return((.self$countIndels() / lenth(.self$cigar_ops))*100)
  },
  
  findDeletions = function(target_start, target_end, alns, del_chars = c("D", "N"),
                           verbose = TRUE ){
    # Get coordinates for narrowing cigars.  For reads with a deletion
    # spanning one or both ends of the target location, narrow the 
    # cigar to encompass the deletion
    #
    # Deletions may be coded as either "D" or "N" (splice junction), 
    # depending upon the mapping software        
      
    # IS LONG DEL ACTUALLY NEEDED?
    
    idxs <- rep(1:length(.self$ref_ranges), lapply(.self$ref_ranges, length))
    genomic <- shift(.self$ref_ranges, start(alns)-1) 
    on_target <- unlist(start(genomic) <= target_end & end(genomic) >= target_start)
    codes <- paste0(idxs, on_target)
    s_del <- !duplicated(codes) & on_target & unlist(.self$cigar_ops) %in% del_chars
    e_del <- rev(!duplicated(rev(codes))) & on_target & unlist(.self$cigar_ops) %in% del_chars
    
    # Get the ranges for narrowing in read coordinates
    result_s <- target_start - (start(alns) - 1)
    result_e <- target_end - target_start  + result_s  
    result_s[idxs[s_del]] <- unlist(start(.self$ref_ranges))[s_del] - 1
    result_e[idxs[e_del]] <- unlist(end(.self$ref_ranges))[e_del] + 1
    
    long_dels <<- s_del | e_del
    return(list(starts = result_s, ends = result_e)) 
  },
  
  getInsertionSeqs = function(){ 
    # Note that the start of a ref_ranges insertion is its genomic end (rightmost base)
    
    ins <- .self$cigar_ops == "I"  
    idxs <- rep(1:length(.self$cigar_ops), sum(ins))
    tseqs <- as.character(mcols(.self$alns)$seq)[idxs]    
    
    if (length(tseqs) == 0) {
        insertions <<- data.frame()
        ins_key <<- vector()
        return()
    }    
    
    qranges <- unlist(.self$query_ranges[ins]) 
    ins_seqs <- as.character(subseq(tseqs, start(qranges), end(qranges)))
    ins_starts <- start(unlist(.self$ref_ranges[ins]))
    genomic_starts <- unlist(start(.self$genome_ranges[ins])) -1 # -1 for leftmost base
    
    df <- data.frame(start = ins_starts, seq = ins_seqs, genomic_start = genomic_starts)
    df$seq <- as.character(df$seq)
    insertions <<- aggregate(rep(1, nrow(df)), by = as.list(df), FUN = table)
    colnames(insertions) <<- c("start", "seq", "genomic_start", "count")
    # Store a key to match the sequences to their insertion
    ins_key <<- match(interaction(df), interaction(insertions[,c(1:3)]))
    names(ins_key) <<- idxs
  },
  
  getInsertionsTable = function(){
    ins_starts <- .self$insertions$genomic_start[.self$ins_key]
    
    if (length(ins_starts) == 0){
        ins <- data.frame()
    } else {
      ins <- data.frame(read = as.integer(names(.self$ins_key)), chromosome = chrom, 
                        start = ins_starts, end = ins_starts, 
                        inserted = .self$insertions$seq[.self$ins_key], stringsAsFactors = FALSE)
    }
    
  },
  
  
  .checkNonempty = function(){
    if (class(.self$alns) == "uninitializedField"){
      message("No on target alignments")
      return(FALSE)
    }
    return(TRUE)
  },
  
  getVariants = function(ref_genome, chrom = NULL, ensembl = FALSE, strand = "+"){
    # Get a data frame of variants and their coordinates for predicting effects,
    # one variant per line
    
    # RefGenome is a BSGenome obj (or FaFile or FaFileList TEST THIS)
    # Return value "read" is the index of the alignment with the variant in .self$alns 
    
    # Strand is wrt the target - "+" for reference strand
    
    if (! .self$.checkNonempty()){
      return(list())
    }
    
    # Get the reference sequence
    if (is.null(chrom)){
      chrom <- as.character(seqnames(.self$alns)[1])
    }
    if (! grepl("chr", chrom)){       
      chrom <- paste0("chr", chrom)
    }
    
    get_start <- min(start(.self$alns))
    get_end <- max(end(.self$alns))
    ref_seq <- getSeq(ref_genome, chrom, get_start, get_end)
    ref_offset <- start(.self$alns) - get_start
    
    if (strand == "-"){  # alignments and ranges have already been reverse complemented
     ref_seq <- reverseComplement(ref_seq)
     ref_offset <- get_end - end(.self$alns) 
    }
        
    # Get deletions:
    dels <- crun$cigar_ops %in% c("D", "N")
    del_idxs <- rep(seq_along(dels), sum(dels))
    
    if (length(del_idxs) == 0){
      dels <- data.frame()
    }else{
      del_ranges <- unlist(.self$ref_ranges[dels]) 
      del_seqs <- as.character(Views(ref_seq, shift(del_ranges, ref_offset[del_idxs])))
      del_gen <- unlist(.self$genome_ranges[dels])
      dels <- data.frame(read = del_idxs, chromosome = chrom, start = start(del_gen), 
                       end = end(del_gen), deleted = del_seqs, stringsAsFactors = FALSE)
    }

    # Get insertions: already computed at initialisation
    ins_starts <- .self$insertions$genomic_start[.self$ins_key]
    
    if (length(ins_starts) == 0){
        ins <- data.frame()
    } else {
      ins <- data.frame(read = as.integer(names(.self$ins_key)), chromosome = chrom, 
                      start = ins_starts, end = ins_starts, 
                      inserted = .self$insertions$seq[.self$ins_key], stringsAsFactors = FALSE)
    }
    
    # Get mismatches - check for mismatches and the read ranges with "M" (match/mismatch)
    matches <- lapply(.self$cigar_ops, function(x) which(x == "M"))
    match_ranges <- .self$query_ranges[matches]
    key <- rep(1:length(match_ranges), lapply(match_ranges, length))
    match_ranges <- unlist(match_ranges)
    match_seqs <- subseq(mcols(.self$alns)$seq[key], start(match_ranges), end(match_ranges))
    
    ref_ranges_ <- unlist(shift(.self$ref_ranges[matches], ref_offset))
    ref_seqs <- as(Views(ref_seq, ref_ranges_), "DNAStringSet")    
    idxs <- which(ref_seqs != match_seqs)  
     
    if (strand == "-"){
      # Make sequences read from genomic left to right for ease of getting genomic coords  
      match_seqs <- DNAStringSet(lapply(match_seqs, rev))
      ref_seqs <- DNAStringSet(lapply(ref_seqs, rev))
    }

    gen_starts <- start(unlist(.self$genome_ranges[matches]))
    
    mm <- do.call(rbind, lapply(idxs, function(i){
            r1 <- as.matrix(ref_seqs[i])
            m1 <- as.matrix(match_seqs[i])
            r1[m1 == "N"] <- "N"
            neq <- m1 != r1 
            if (length(which(neq)) == 0){ # Mismatches are only due to "N"
              return(data.frame())
            }
            gen_starts_i <- gen_starts[i] + (which(neq) -1 )
            df <- data.frame(read = key[i], chromosome = chrom, 
                             start =gen_starts_i, ref = r1[neq], query = m1[neq])
            
          }))
    
    vars <- list(insertions = ins, deletions = dels, mismatches = mm)  
    
    if (ensembl == TRUE){
        return(getVarsEnsemblFormat(vars, strand))
    } else{
      return(vars)
    }
  },

  getVarsEnsemblFormat = function(vars, strand = "*"){
    # Ensembl requires end = start - 1 for insertions
    # End points are included
    
    if (! .self$.checkNonempty()){
      return(list())
    }
        
    vins <- vars$insertions
    ins <- data.frame(chrom = vins$chromosome, start = vins$start, end = vins$end -1 , 
             allele = sprintf("-/%s", vins$inserted), id = vins$read)
    
    vdels <- vars$deletions
    dels <- data.frame(chrom = vdels$chromosome, start = vdels$start, end = vdels$end,
                       allele = sprintf("%s/-", vdels$deleted), id = vdels$read)
    
    vmm <- vars$mismatches
    mm <- data.frame(chrom = vmm$chromosome, start = vmm$start, end = vmm$start,
                     allele = sprintf("%s/%s", vmm$ref, vmm$query), id = vmm$read)
                     
    vars <- rbind(ins,dels,mm)
    if (nrow(vars) > 0){
      vars$strand <- strand
      vars <- vars[order(vars$id),]
      vars$id <- names(.self$alns[vars$id])
      vars <- vars[,c("chrom","start","end","allele","strand","id")]
    }
    vars$chrom <- gsub("chr", "", vars$chrom)
    vars
  },
  
  getCigarLabels = function(short = TRUE, match_label = "no variant", 
                            genome_to_target = NULL, rc = FALSE, target_start = NULL,
                            target_end = NULL, ref = NULL, split_non_indel = TRUE, 
                            mismatch_label = "SNV", cut_site = 18){
    
    # Sets / returns cigar labels                        
    # If short = only consider insertions and deletions (except for wildtype?)
                            
    # genome_to_target, if provided, is a vector with names being genomic start coords
    # and values being coords with respect to the target site  
    
    # Shorten? -> Renumber?  If not renumbered, record starting location when different from 
    # target start
    
    if (! class(.self$cigar_labels) == "uninitializedField"){
      return(.self$cigar_labels)
    }
    
    # for default cigar (short = FALSE) and short cigar without renumbering:  
    # indicate the starting location with respect to the target start if it is not the target
    start_offset <- rep("", length(.self$alns))
    
    if (short == FALSE | is.null(genome_to_target)){
      if (rc == TRUE){
        # Note: genome_ranges for rc are written rightmost to leftmost
        
        start_offset <- target_end - unlist(lapply(.self$genome_ranges, function(x) max(end(x))),
                                            use.names = FALSE)
        start_offset_new <- target_end - max(end(.self$genome_ranges))
        #print(c(start_offset_new, start_offset, start_offset_new == start_offset))
        
      } else {
        start_offset <- unlist(lapply(.self$genome_ranges, function(x) min(start(x))), 
                               use.names = FALSE) - target_start
        start_offset_new <-  min(start(.self$genome_ranges)) - target_start 
        #print(start_offset_new == start_offset)                      
      }
      is.nonzero <- start_offset != 0
      start_offset[is.nonzero] <- sprintf("(%s)", start_offset[is.nonzero])
      start_offset[! is.nonzero] <- ""       
      
      if ( short == FALSE ){
        cigar_labels <<- paste0(start_offset, cigar(.self$alns))
        return(.self$cigar_labels)
      }
    }
    
    #idxs <- .self$cigar_ops != "M"
    #ops <- .self$cigar_ops[idxs]
    #rranges <- .self$ref_ranges[idxs]
    #qranges <- .self$query_ranges[idxs]
    #cig_idxs <- rep(1:length(idxs), lapply(idxs, length))
    #cig_idxs <- .self$long_dels[rep(1:length(idxs), sum(idxs == TRUE))]
       
    idxs <- lapply(.self$cigar_ops, function(x) which(x != "M"))
    ops <- lapply(seq_along(idxs), function(i) .self$cigar_ops[[i]][idxs[[i]]])
    rranges <- .self$ref_ranges[idxs]
    qranges <- .self$query_ranges[idxs]
    cig_idxs <- rep(1:length(idxs), lapply(idxs, length))   
       
         
    if (is.null(genome_to_target)){
      start <- unlist(start(rranges))
      end <- unlist(end(rranges))
    } else {
      gen_ranges <- .self$genome_ranges[idxs]     
      if (rc == FALSE){
        start <- genome_to_target[as.character(unlist(start(gen_ranges)))]
        end <- genome_to_target[as.character(unlist(end(gen_ranges)))]
    
      }else {
        start <- genome_to_target[as.character(unlist(end(gen_ranges)))]
        end <- genome_to_target[as.character(unlist(start(gen_ranges)))]
      }
    }
     
    info <- data.frame(op = unlist(ops), qwidth = unlist(width(qranges)), 
              start = start, end = end, rwidth = unlist(width(rranges)), 
              long_del = .self$long_dels[cig_idxs])
              
    getShortOp <- function(mutn){
       if (mutn["op"] %in% c("N", "D")){ 
         sprintf('%s:%sD', as.numeric(mutn["start"]), as.numeric(mutn["rwidth"]))
       }else if (mutn["op"] == "I"){
         sprintf('%s:%sI', as.numeric(mutn["start"]), as.numeric(mutn["qwidth"]))
       }
       
    }
    result <- rep(match_label, length(idxs))
    if (nrow(info) == 0) return(result)
    short_ops <- apply(info, 1, getShortOp)
    short_ops <- split(short_ops, cig_idxs)
    short_ops <- sapply(short_ops, function(x) do.call(paste, c(as.list(x),sep = ",")))
    result[as.numeric(names(short_ops))] <- short_ops
    if (split_non_indel){
      result <- .splitNonIndel(ref, result, match_label, mismatch_label, cut_site)
    }
    cigar_labels <<- result
    result
  },
  
  .splitNonIndel = function(ref, cig_labels, match_label = "no variant", 
                             mismatch_label = "SNV", cut_site = 18){
    no_var <- which(cig_labels == match_label & mcols(.self$alns)$seq != ref)
    if (length(no_var) == 0) return(cig_labels)
        
    no_var_seqs <- as.matrix(mcols(.self$alns)$seq[no_var])
    rr <- strsplit(as.character(ref), "")[[1]]
    result <- apply(no_var_seqs, 1, function(x) sprintf("%s:%s", mismatch_label, 
                     paste(which((x != rr & x != "N")) - cut_site -1, collapse = ",")))
    result[result == sprintf("%s:", mismatch_label)] <- match_label
    cig_labels[no_var] <- result 
    return(cig_labels)       
  }
  
)

#_______________________________________________________________________________________


CrisprSet = setRefClass(
  Class = "CrisprSet",
  fields = c("crispr_runs", 
             "ref",
             "insertion_sites",
             "nonempty_runs",
             "cigar_freqs",
             "target",
             "genome_to_target",
             "pars")
)

CrisprSet$methods(
  initialize = function(bam_fnames, ref, target_chr, target_start, target_end, rc = FALSE,
                        short_cigars = FALSE, names = NULL, exclude_ranges = GRanges(), 
                        exclude_names = NA, merge_chimeras = TRUE, renumbered = FALSE, 
                        target_loc = NA, match_label = "no variant", verbose = TRUE, ...){
    
    # TO DO - DECIDE WHETHER TO KEEP EMPTY RUNS
    strand <- ifelse(rc == TRUE, "-", "+")
    target <<- GRanges(target_chr, IRanges(target_start, target_end), strand = strand)    
    ref <<- ref 
    pars <<- list("match_label" = match_label, "target_loc" = target_loc, 
                  "mismatch_label" = "SNV", "renumbered" = renumbered)
    
    #pars <<- modifyList(pars, ...)   
        
    if (verbose == TRUE) cat("Reading bam files\n")
    crispr_runs <<- lapply(bam_fnames, CrisprRun$new, target_chr, target_start, 
                           target_end, rc = rc, exclude_ranges = exclude_ranges,
                           exclude_names = exclude_names, merge_chimeras = merge_chimeras,
                           verbose = verbose, ...) 
    
    if (! is.null(names)) names(.self$crispr_runs) <- names
    nonempty_runs <<- sapply(.self$crispr_runs, function(x) {
                             ! class(x$alns) == "uninitializedField"})
    
    .self$crispr_runs <<- .self$crispr_runs[.self$nonempty_runs]
    if (length(.self$crispr_runs) == 0) stop("no on target runs")
    
    if (verbose == TRUE) cat("Renaming cigar strings and counting indels\n")
    cig_by_run <- .self$.setCigarLabels(renumbered = renumbered, target_loc = target_loc,
                           target_start = target_start, target_end = target_end, 
                           rc = rc, match_label = match_label, ref = ref)
    .self$.countCigars(cig_by_run)
  },
  
  show = function(){
    print(c(class(.self), sprintf("CrisprSet object containing %s CrisprRun samples", 
            length(.self$crispr_runs)), .self$cigar_freqs))
  },
  
  .setCigarLabels = function(renumbered = FALSE, target_loc = NA, target_start = NA,
                            target_end = NA, rc = FALSE, match_label = "no variant",
                            ref = NULL){
      g_to_t <- NULL
      
      if (renumbered == TRUE){
        if (any(is.na(c(target_loc, target_start, target_end, rc)))){
          stop("Must specify target_loc (cut site), target_start, target_end and rc
                for renumbering")
        }
        g_to_t <- genomeToTargetLocs(target_loc, target_start, target_end, rc)
      }
      
      cut_site <- ifelse(is.na(target_loc), 18, target_loc)
      
      cig_by_run <- lapply(.self$crispr_runs,
            function(crun) crun$getCigarLabels(match_label = match_label, rc = rc,
                                               genome_to_target = g_to_t, ref = ref, 
                                               cut_site = cut_site))     
      return(cig_by_run)
  }, 
   
  .countCigars = function(cig_by_run = NULL, nonvar_first = TRUE){
    # Note that this function does not consider starts, two alignments starting at
    # different locations but sharing a cigar string are considered equal
    
    if (is.null(cig_by_run)){
      cig_by_run <- lapply(.self$crispr_runs, function(crun) crun$getCigarLabels())
    }
    
    unique_cigars <- unique(unlist(cig_by_run))  
    
    m <- matrix(unlist(lapply(cig_by_run, function(x) table(x)[unique_cigars])), 
          nrow = length(unique_cigars), dimnames = list(unique_cigars, names(.self$crispr_runs)))
    
    m[is.na(m)] <- 0
    m <- m[order(rowSums(m), decreasing = TRUE),, drop = FALSE]
    
    if (nonvar_first){    
      is_ref <- grep(.self$pars["match_label"], rownames(m))
      is_snv <- grep(.self$pars["mismatch_label"], rownames(m))
      new_order <- setdiff(1:nrow(m), c(is_ref,is_snv))
      m <- m[c(is_ref, is_snv, setdiff(1:nrow(m), c(is_ref,is_snv))),,drop = FALSE]
    }
    
    cigar_freqs <<- m
  },
  
  .getFilteredCigarTable = function(top_n, freq_cutoff){
    rs <- rowSums(.self$cigar_freqs)
    minfreq <- rs >= freq_cutoff
    topn <- rank(-rs) <= top_n
    cig_freqs <- .self$cigar_freqs[minfreq & topn ,, drop = FALSE] 
    return(cig_freqs) 
  },
  
  countVariantAlleles = function(counts_t = NULL){
    if (is.null(counts_t)) counts_t <- .self$cigar_freqs
    counts_t <- counts_t[!rownames(counts_t) == .self$pars["match_label"],,drop = FALSE]
    alleles <- colSums(counts_t != 0)    
    return(data.frame(Allele = alleles, Sample = names(alleles)))
  },
  
  heatmapCigarFreqs = function(as_percent = FALSE, x_size = 16, y_size = 16, 
                               x_axis_title = NULL, x_angle = 90, annotate_counts = TRUE, 
                               freq_cutoff = 0, top_n = nrow(.self$cigar_freqs), ...){
    
    cig_freqs <- .getFilteredCigarTable(top_n, freq_cutoff)
    p <- cigarFrequencyHeatmap(cig_freqs, as_percent, x_size, y_size, x_axis_title,
                               x_angle, annotate_counts, ...)
    return(p)
  },
  
  plotVariants = function(freq_cutoff = 0, top_n = nrow(.self$cigar_freqs), 
                          short_cigars = FALSE, renumbered = .self$pars["renumbered"], 
                          show_genomic = FALSE, ...){
                          
    # var_freq_cutoff = i (integer) only plot variants that occur >= i times
    # top_n = total number of variants to plot
    # ... arguments for plotAlignments
    
    cig_freqs <- .getFilteredCigarTable(top_n, freq_cutoff)
  
    alns <- .self$makePairwiseAlns(cig_freqs)
    if (class(.self$insertion_sites) == "uninitializedField" | 
        !("cigar" %in% colnames(.self$insertion_sites))){
      .self$getInsertions() 
    }
    
    # How should the x-axis be numbered? 
    # Baseline should be numbers, w optional genomic locations
    if (renumbered == TRUE){
      genomic_coords <- c(start(.self$target):end(.self$target))
      target_coords <- .self$genome_to_target[as.character(genomic_coords)]
      if (as.character(strand(.self$target)) == "-"){
        target_coords <- rev(target_coords)
      }
      xbreaks = which(target_coords %% 5 == 0 | abs(target_coords) == 1)
      target_coords <- target_coords[xbreaks]
      
      p <- plotAlignments(.self$ref, alns, .self$insertion_sites, 
                           xtick_labs = target_coords, xtick_breaks = xbreaks, ...)
    } else {
      p <- plotAlignments(.self$ref, alns, .self$insertion_sites, ...)    
    }
    return(p)
  },
  
  getInsertions = function(with_cigars = TRUE){
    if (with_cigars == FALSE){
      all_ins <- do.call(rbind, lapply(.self$crispr_runs, function(x) x$insertions))
    } else {
      all_ins <- do.call(rbind, lapply(.self$crispr_runs, function(x) {
                ik <- x$ins_key
                v <- data.frame(ik, x$getCigarLabels()[as.integer(names(ik))])
                v <- v[!duplicated(v),]
                v <- v[order(v$ik),]
                cbind(x$insertions[v[,1],], cigar = v[,2])
              }))
    }
    if (nrow(all_ins) == 0) {
      insertion_sites <<- all_ins
      return()
    }
    insertion_sites <<- all_ins[order(all_ins$start, all_ins$seq),, drop = FALSE]
  },
  
  makePairwiseAlns = function(cig_freqs = .self$cigar_freqs, ...){
    # Get alignments by cigar string, make the alignment for the consensus
    # The short cigars (not renumbered) do not have enough information, 
    # use the full cigars for sorting

    cigs <- unlist(lapply(.self$crispr_runs, function(x) cigar(x$alns)), use.names = FALSE)
    cig_labels <- unlist(lapply(.self$crispr_runs, function(x) x$getCigarLabels()), use.names = FALSE)
    
    names(cigs) <- cig_labels # calling by name with duplicates returns the first match
      
    splits <- split(seq_along(cig_labels), cig_labels)
    splits <- splits[match(rownames(cig_freqs), names(splits))]
    
    splits_labels <- names(splits)
    names(splits) <- cigs[names(splits)]
        
    x <- lapply(.self$crispr_runs, function(x) x$alns)
    all_alns <- do.call(c, unlist(x, use.names = FALSE))
  
    seqs <- c()
    starts <- c()
    
    # SOMEWHERE HERE - CONSENSUS SEQ ONLY TAKING ONE SAMPLE?
        
    for (i in seq_along(splits)){
      idxs <- splits[[i]]
      seqs[i] <- consensusString(mcols(all_alns[idxs])$seq)
      start <- unique(start(all_alns[idxs]))
      if (length(start) > 1)
        stop("Sequences with the same cigar string have different starting locations.
              This case is not implemented yet.")  
      starts[i] <- start[1]
    }  
   
   alns <- mapply(seqsToAln, names(splits), seqs, aln_start = starts, 
                target_start = start(.self$target), target_end = end(.self$target), ...)
   
   names(alns) <- splits_labels
   alns
  },
  
  plotVariantOverview = function(){
  
    # heatmap must take same filtering args, OR be able to match the alignment names
  
  },
  
  genomeToTargetLocs = function(target_loc, target_start, target_end, rc = FALSE){
    # target_loc should be relative to the start of the target sequence, even if the 
    # target is on the negative strand
    # target_loc is the left side of the cut site (Will be numbered -1)
    # target_start and target_end are genomic coordinates, with target_start < target_end
    # rc: is the target on the negative strand wrt the reference?
    # returns a vector of genomic locations and target locations
    
    # Example:  target_loc = 5
    # Before: 1  2  3  4  5  6  7  8 
    # After: -5 -4 -3 -2 -1  1  2  3
    # Left =  original - target_loc - 1
    # Right = original - target_loc
        
    gs <- min(unlist(lapply(.self$crispr_runs, function(x) start(unlist(x$genome_ranges)))))
       
    all_gen_ranges <- lapply(.self$crispr_runs, function(x) x$genome_ranges)
    
    # Nope - here gives a vector
    #gs_test <- min(start(do.call(c, unlist(all_gen_ranges, use.names = FALSE))))
    #cat(sprintf("testing gs: %s = %s: %s", gs_test, gs, gs_test == unlist(gs, use.names=FALSE)))
    
    ge <- max(unlist(lapply(.self$crispr_runs, function(x) end(unlist(x$genome_ranges)))))
    
    if (rc == TRUE){
      tg <- target_end - (target_loc - 1)
      new_numbering <- rev(c(seq(-1*(ge - (tg -1)),-1), c(1:(tg - gs))))
      names(new_numbering) <- c(gs:ge)
      
    } else {
      tg <- target_start + target_loc - 1
      new_numbering <- c(seq(-1*(tg - (gs-1)),-1), c(1:(ge - tg)))
      names(new_numbering) <- c(gs:ge)
    }
    genome_to_target <<- new_numbering
    new_numbering
  }
)


#__________________________________________________________


CrisprMultiplex = setRefClass(
  Class = "CrisprMultiplex"#,
  #fields = c()
)


CrisprMultiplex$methods(

  initialize = function(bam_fnames, crispr_targets, names = NA, pcr_targets = NA, 
                        merge_chimeras = FALSE, exclude_chimeras = TRUE, 
                        exclude_ranges = GRanges(), exclude_names = NA, mapq = NA){
                        
    alns_by_tg_by_bam <- mapply(readMultiplexBam, bam_fnames, crispr_targets)  
    
     
  },
    
  readMultiplexBam = function(bam_fname, crispr_targets, pcr_targets = NULL, tolerance = 4,
                              merge_chimeras = FALSE, exclude_chimeras = TRUE, 
                              tag_chimeras = FALSE, name = NA, exclude_ranges = GRanges(),  
                              exclude_names = NA, mapq = NA, verbose = TRUE){
    
    # mapq - remove reads with quality equal or below this value
    # allow pcr_tolerance on either side of the pcr target
    # targets: GRanges 
    # ReadGAlignments doesn't return unmapped reads, isUnmappedQuery doesn't work
    # regions / alignments to exclude are excluded prior to detecting chimeras
    # (a chimera consisting of two pieces, one in an excluded region, 
    # is therefore not counted as chimeric)
    
    if (is.na(mapq)){
      param <- ScanBamParam(what = c("seq"))
      bam <- readGAlignments(bam_fname, param = param, use.names = TRUE)
      if (verbose == TRUE) ninitial <- length(bam)
    } else {
      param <- ScanBamParam(what = c("seq", "mapq"))
      bam <- readGAlignments(bam_fname, param = param, use.names = TRUE)
      if (verbose == TRUE) ninitial <- length(bam)
      bam <- bam[mcols(bam)$mapq > mapq]
    }
    
    if (verbose == TRUE){
      cat(sprintf("Read %s alignments %s\n", ninitial, ifelse(is.na(mapq), "",
      sprintf("(%s with mapq < %s excluded, %s remaining)", ninitial - length(bam), 
      mapq, length(bam)))))
    }
    
    # TO DO HERE: EXCLUDE NAMES, EXCLUDE REGIONS    
    
    
    
    if (any(c(merge_chimeras,exclude_chimeras,tag_chimeras)) == TRUE){
      bam <- classifyChimeras(bam, exclude = exclude_chimeras, merge = merge_chimeras, 
                            tag = tag_chimeras, verbose = verbose, name = name)
    }
    
    if (verbose == TRUE) cat("%s reads after filtering chimeras\n", length(bam))
    
    
    # Distribute reads between targets
    
    if (! is.null(pcr_targets)){        
    
      # Extrapolate genomic coords for clipped sections
      l_clip <- as.numeric(gsub(".*[A-Z].*", 0, gsub("[HS].*","", cigar(bam))))
      r_clip <- as.numeric(gsub("^$", 0, gsub('.*M|[HS]$', "", cigar(bam)))) 
      bamgr <- GRanges(seqnames(bam), IRanges(start(bam) - l_clip, end(bam) + r_clip))
      
      hits_pcr <- findOverlaps(bamgr, pcr_targets, ignore.strand = TRUE, 
                               type = "equal", maxgap = tolerance)
      
      hits_pcr_l <- length(unique(hits_pcr@queryHits))
      if (verbose == TRUE){
        cat(sprintf("%s from %s (%.2f%%) reads overlap a pcr region\n", 
            hits_pcr_l, length(bam), hits_pcr_l/length(bam)*100))
        remaining <- setdiff(c(1:length(bam)), hits_pcr@queryHits) 
        rhits <- findOverlaps(bamgr[remaining], pcr_targets)  
        rhitsl <- length(unique(rhits@queryHits))   
        cat(sprintf("Of the remaining %s reads, %s (%.2f%%) partially overlap a pcr region\n\n", 
            (length(bam)-hits_pcr_l),rhitsl, rhitsl/(length(bam)-hits_pcr_l)*100))
      }
      
      if (any(duplicated(hits_pcr@queryHits))){
        warning("Cannot distinguish pcr targets with this tolerance")
      }      
    }
    
    
    ## Find reads that partially cover the target site 
    #partially_on_target <- factor(as.character((start(bam) <= target_start | end(bam) >= target_end) & 
    #                 seqnames(bam) == target_chr), levels = c("TRUE", "FALSE"))
    #partial_names <- split(names(bam), partially_on_target)    
    #names(partial_names) <- c("partial", "off_target")  
       
    ## Split reads into fully on-target and off-target    
    #is_on_target <- factor(as.character(start(bam) <= target_start & end(bam) >= target_end & 
    #                 seqnames(bam) == target_chr), levels = c("TRUE", "FALSE"))


    by_target <- split(queryHits(rhits), subjectHits(rhits))

  }   
  
)
