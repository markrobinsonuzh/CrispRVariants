library(GenomicAlignments)
library(Rsamtools)
library(ggplot2)
library(reshape2)
library(gridExtra)

# sangerseqR, GenomeFeatures should be in "sugests"

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
# To do - deal with non-ATGC characters
# To do - indel realignment around target loc
# To do - make read classification optional?
# Can I assume that all reads in a set have different names?
# Warning that normal cigar misleading in aln plots
# Take is.na out of plotAlignments?
# automatic layout heights and widths - absolute for gene track?


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


writeFastq <- function(outf, vals){
    o <- file(outf, "a")
    seqname <- sprintf("@%s", vals$seqname)
    qualname <- sprintf("+%s", vals$seqname)
    writeLines(c(seqname, vals$seq, qualname, vals$quals) , o)
    close(o)
}

abifToTrimmedFastq <- function(seqname, fname, outfname, trim = TRUE, cutoff = 0.05, 
                               min_seq_len = 20, offset = 33){
        
    sangerseqr <- require(sangerseqR)
    stopifnot(sangerseqr == TRUE)
    
    # Translation of python function
    abif <- read.abif(fname)
    if (is.null(abif@data$PCON.2)){
       print(sprintf("failed on %s", seqname))
       return()
    }
    
    # Remove the extra character if it exists
    nucseq <- substring(abif@data$PBAS.2, 1, length(abif@data$PLOC.2))
    num_quals <- utf8ToInt(abif@data$PCON.2)[1:length(abif@data$PLOC.2)] 
     
    if (trim == FALSE){
        writeFastq(outfname, list("seq" = nucseq, "quals" = rawToChar(as.raw(num_quals + 33))))
        return()
    }
    
    if (nchar(nucseq) <= min_seq_len){
        stop('Sequence can not be trimmed because it is shorter than the trim segment size')
    } 
    scores = cutoff - 10^(num_quals / -10)
    running_sum <- rep(0, length(scores) + 1)
     
    # Note running_sum counts from zero, scores from 1 
    for (i in 1:length(scores)){
       num <- scores[i] + running_sum[i]
       running_sum[i+1] <- ifelse(num < 0, 0, num)
    }

    trim_start <- min(which(running_sum > 0)) - 1
    trim_finish <- which.max(running_sum) - 2 
    # -1 for running_sum offset, -1 because python doesn't include ends

    writeFastq(outfname, list("seqname" = seqname, "seq" = substring(nucseq, trim_start, trim_finish),
                "quals" = rawToChar(as.raw(num_quals[trim_start:trim_finish]+offset))))
    return()
}

reverseCigar <- function(cigar){
   cigar.widths <- rev(strsplit(cigar, '[A-Z]')[[1]])
   cigar.ops <- rev(explodeCigarOps(cigar)[[1]])
   paste0(cigar.widths,cigar.ops, collapse = "")
}

seqsToAln <- function(cigar, dnaseq, del_char = "-", aln_start = NULL, target_start = NULL, 
                      target_end = NULL){
    # Remove insertion sequences
    # When genomic coordinates for the dnaseq alignment start and the target region are
    # provided, aligned sequence is cropped to the target region if necessary
        
    wrt_ref <- cigarRangesAlongReferenceSpace(cigar)[[1]]
    wrt_read <- cigarRangesAlongQuerySpace(cigar)[[1]]
    ops <- explodeCigarOps(cigar)[[1]]
    segs <- as.character(Views(dnaseq, wrt_read))
    segs[which(ops == "I")] <- ""
    for (j in which(ops == "D")){
         segs[j] <- paste0(rep(del_char, width(wrt_ref[j])), collapse = "")
      }
    result <- paste0(segs, collapse = "")
    
    if (! is.null(aln_start) & ! is.null(target_start) & ! is.null(target_end) ){
      trim_start <- target_start - (aln_start - 1)

      if ( trim_start < 0){
        stop("dnaseq to be trimmed must start before the target location")
      }
            
      trim_end <- trim_start + (target_end - target_start)
      if (trim_end < length(result)){
        stop("dnaseq is not long enough to trim to the target region")
      }
      result <- subseq(result, trim_start,trim_end)
    }
    result
}

cigarFrequencyHeatmap <- function(cigar_freqs, as_percent = TRUE, x_size = 10, 
                             x_axis_title = NULL, x_angle = 90, annotate_counts = TRUE,
                             col_sums = TRUE, colours = "yellowred"){
  # add_blank - should a blank row be added to the top?  
  #(useful for aligning heatmap with variants plot) if col_sums = FALSE
    
  cig_freqs <- cigar_freqs
  if (col_sums == TRUE){
    # Make space for totals to be added
    cig_freqs <- rbind(Total = rep(NA, ncol(cigar_freqs)), cigar_freqs)
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
    g <- g + geom_text(data = counts, aes(label = Count, fill = NULL, fontface = ff))
  }
  g <- g + ylab(NULL) + xlab(x_axis_title) + theme_bw() +
       theme(axis.text.x = element_text(size = x_size, angle = x_angle))
        
  if (colours == "yellowred"){     
    hmcols<-colorRampPalette(c("white","gold","orange","orangered","red", "darkred"))(50) 
    g <- g + scale_fill_gradientn(colours = hmcols, na.value = "white") 
  }            
  return(g)
}


getCodonFrame <- function(txdb, target_chr, target_start, target_end){
  # Target_chr must match txdb
  require(VariantAnnotation)
  refLocsToLocalLocs(GRanges(target_chr, IRanges(target_start, target_end)), txdb)
  
}

nucleotideToAA <- function(seqs, txdb, target_start, target_end){
    
}

plotAlignments <- function(ref, alns, ins_sites, pam_loc = NA, show_plot = FALSE, target_loc = 18,
                           pam_start = NA, pam_end = NA, ins_size = 4, legend_cols = 3,
                           xlab = NULL){
  
  # ref: the reference sequence
  # alns: named vector of aligned sequences, with insertions removed
  # ins_sites:  table of insertion_sites, must include cols named "start", "cigar" and "seq"
  # pam_loc: location of PAM with respect to the target site
  # Insertion locations are determined by matching ins_sites$cigar with names(alns)
  # All characters other than ACTG are labelled N
    
  # Reverse alignment order, as ggplot geom_tile plots bottom up
  aln_chrs <- strsplit(c(rev(alns), Reference = as.character(ref)), "")
  temp <- t(as.data.frame(aln_chrs))
  rownames(temp) <- names(aln_chrs)
  m <- melt(temp)
  ambig <- which(! m$value %in% c("A", "C", "T", "G", "N", "-"))
  m$value <- as.character(m$value)
  m[ambig, "value"] <- "N"  
  m$value <- factor(m$value, levels = c("A", "C", "T", "G", "N", "-"))  
  m$isref <- as.character(ifelse(m$Var1 == "Reference", 1, 0.75))
  m_cols <- c("#e41a1c", "#377eb8", "#4daf4a", "#000000", "#CCCCCC","#FFFFFF", "#FFFFFF")
  names(m_cols) <- c("A", "C", "T", "G", "N","-", "+")
  m$cols <- m_cols[m$value]

  # Colours and shapes for the insertion markers
  shps <- c(21,23,25) 
  clrs <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
            "#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#661100",
            "#CC6677","#882255", "#AA4499")
                                        
  # Make a data frame of insertion locations, format seqs to fixed width for nicer legends
  ins_ord <- match(ins_sites$cigar, names(aln_chrs))  
  seqs <- ins_sites[!is.na(ins_ord),"seq"]
  max_seq_ln <- max(sapply(seqs, nchar)) + 2 
  seqs <- sprintf(paste0("%-",max_seq_ln,"s"), seqs)
  
  ins_points <- data.frame(x = ins_sites[!is.na(ins_ord),"start"] - 0.5,
                           y = na.omit(ins_ord) + 0.45,
                           seq = seqs)

  ins_points$shapes <- as.factor(1:nrow(ins_points))
  ins_points$colours <- as.factor(rep(clrs, 3)[1:nrow(ins_points)])
  fill_clrs <- rep(clrs, 3)[1:nrow(ins_points)]
  fill_shps <- rep(shps, 17)[1:nrow(ins_points)]
  legend_nrow <- ceiling(nrow(ins_points)/ legend_cols)

  tile_height <- 0.5
  p <- ggplot(m, aes(x = Var2, y = Var1, fill = cols))+
     geom_tile(aes(alpha = isref), height = tile_height)+ 
     geom_text(aes(label = value), size = 4) + 
     scale_alpha_manual(values = c(0.5,1), guide = "none") + 
     #scale_x_continuous(expand = c(0,0.5), breaks = set_breaks, labels = set_vals) + # expand is the distance from the axis, multiplicative + additive
     ylab(NULL) + xlab(xlab) +
     theme_bw() + theme(axis.text.y = element_text(size = 12),axis.text.x = element_text(size = 8),
                        axis.title.x = element_text(vjust = -0.5), legend.position = "bottom")
     
  p <- p + geom_point(data = ins_points, aes(x = x, y = y, shape = shapes, fill = colours),
                      colour = "#000000", size = ins_size)  +
     scale_fill_identity() +
     scale_shape_manual(name = "", values = fill_shps, breaks = ins_points$shapes, labels = ins_points$seq) +
     guides(shape = guide_legend(ncol = legend_cols, override.aes = list(fill = fill_clrs, size = 6))) 
     p <- p + theme(legend.key = element_blank(), legend.text = element_text(size = 16),
                    legend.key.height = unit((legend_nrow * 0.25), "lines"),
                    legend.margin = unit(4, "lines"))
      
  p <- p + geom_vline(xintercept= target_loc + 0.5, colour = "red", linetype = "dotted")
  
  
  
  # If pam_loc is given, highlight the pam in the reference
  #  - 0.5 for tile boundaries not centres
  
  if (! is.na(pam_start)){    
    if (is.na(pam_end)){
      pam_end <- pam_start + 2
    }
    pam_df <- data.frame(xmin = pam_start - 0.5, xmax = pam_end - 0.4,
                         ymin = length(names(aln_chrs)) - (tile_height / 2 + 0.1),
                         ymax = length(names(aln_chrs)) + (tile_height / 2 + 0.1))
    p <- p + geom_rect(data=pam_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax,
                        ymax = ymax, x = NULL, y = NULL),
                         color = "red", size = 1, fill = "transparent") 
   
  }
  
  if (show_plot == TRUE){
    print(p)
  }
  return(p)
}

.tableInsertionSeqs <- function(insertion_points, ncol = 4){
  # insertion_points should have shape, colour, and seq columns
  
  
       #scale_shape_manual(values = ins_points$shapes, labels = ins_points$seq, name = "") +
       #scale_fill_manual(labels = ins_points$seq, values = ins_points$colours, name = "")


}

panelPlot <- function(txdb, target_chr, target_start, target_end, aln_p, heat_p, 
                      col.widths = c(2, 1)){

  
  print("starting panel plot")
  # Make the gene plot
  require(ggbio)
  genes <- genes(txdb)
  target <- GRanges(target_chr, IRanges(target_start, target_end))
  wh <- genes[findOverlaps(genes, target)@queryHits]  
  p1 <- autoplot(txdb, wh)
  
  #Pull off the y limits from the transcript plot
  yranges <- ggplot_build(p1)$panel$ranges[[1]]$y.range  
  xranges <- ggplot_build(p1)$panel$ranges[[1]]$x.range
  target_df <- data.frame(xmin = start(target), xmax = end(target), 
                          ymin = yranges[1], ymax = yranges[2])
    
  p1 <-  p1 + geom_rect(data = target_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                 colour = "red", fill = NA) + theme_bw()

  p1 <- ggplotGrob(as.list(attributes(p1))$ggplot)
  
  heat_p <- heat_p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

  p2 <- ggplotGrob(aln_p)
  p3 <- ggplotGrob(heat_p)

  maxHeight = grid::unit.pmax(p2$heights[2:5], p3$heights[2:5])
  p2$heights[2:5] <- as.list(maxHeight)
  p3$heights[2:5] <- as.list(maxHeight)

  #return(arrangeGrob(p2, p3, ncol = 2, widths = col.widths))

  return(grid.arrange(p1, arrangeGrob(p2, p3, ncol = 2, widths = col.widths), 
         nrow = 2, heights = c(2, 8)))

}





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
             "read_types")
)

CrisprRun$methods(
  initialize = function(bam_fname, target_chr, target_start, target_end, rc = FALSE, 
                        name = NULL, count_unmapped = TRUE, merge_chimeras = FALSE,
                        exclude_ranges = GRanges(), exclude_names = NA, ...){
                        
    # ... args for merging chimeras
    
    name <<- ifelse(is.null(name), bam_fname, name)
    galns <- readBam(bam_fname, target_chr, target_start, target_end, merge_chimeras, 
                     exclude_ranges, exclude_names, ...)
    
    if (count_unmapped == TRUE){
      param <- ScanBamParam(what = c("qname"), flag=scanBamFlag(isUnmappedQuery=TRUE))
      read_types[scanBam(bam_fname, param = param)[[1]]$qname] <<- "unmapped"
    }
        
    off_target <<- galns[["off_target"]]
    on_target <- galns[["on_target"]]
        
    if (length(on_target) == 0) return()

    ref_ranges <<- cigarRangesAlongReferenceSpace(cigar(on_target))
    query_ranges <<- cigarRangesAlongQuerySpace(cigar(on_target))
    
    cigar_ops <<- explodeCigarOps(cigar(on_target)) 
    .self$readsToTarget(on_target, target_start, target_end, rc)
    .self$getInsertionSeqs()
  },
  
  readBam = function(bam_fname, target_chr, target_start, target_end, merge_chimeras = FALSE,
                     exclude_ranges = GRanges(), exclude_names = NA, ...){
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
  
    param <- ScanBamParam(what = c("seq"))
    bam <- readGAlignments(bam_fname, param = param, use.names = TRUE)
    read_types <<- sapply(unique(names(bam)), function(x) NA)
    
    # Exclude by range
    fo <- findOverlaps(bam, exclude_ranges)@queryHits
    read_types[unique(names(bam[fo]))] <<- "excluded"
    bam <- bam[setdiff(seq_along(bam), fo)]
    
    if (length(bam) == 0) return(list("on_target" = NULL, "off_target" = NULL))
    
    # Exclude by name   
    excluden <- which(names(bam) %in% exclude_names)     
    read_types[unique(excluden)] <<- "excluded"
    bam <- bam[setdiff(seq_along(bam), excluden)]     
    
    if (merge_chimeras == TRUE){
      # Split the bam file by seq name, preserving the original order
      bam_by_name <-split(bam, factor(names(bam), levels = unique(names(bam))))
      result <- unlist(lapply(bam_by_name, .self$annotateChimericAlns, ...), use.names = FALSE)
      bam <- bam[result]
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
    
    # Note that a read can currently be classified as "excluded" if it's a chimera
    # and half is excluded
    
    remaining <- names(which(is.na(read_types) | read_types == "excluded"))
    read_types[remaining[remaining %in% names(result$off_target)]] <<- "off_target"
    read_types[remaining[remaining %in% partial_names$partial]] <<- "partially_on_target"
    read_types[remaining[remaining %in% names(result$on_target)]] <<- "on_target"
    result
  },   
  
  annotateChimericAlns = function(aln_set, keep_multichr = TRUE){  
    # Case: Aligned regions overlap
    #  1-2-3-4-5
    #      3-4-5-6-7
    #
    # Case: Inversion:
    # 1-2-3-4-5
    #            6-5-4-3 
    
    # keep_multichr: should members of a multi-chromosome chimera be excluded?
    # note that only the on-target members will be counted
    
    
    # Chimeras must have multiple seqs on the same chromosome and strand
    result <- rep(TRUE, length(aln_set))
    
    if (length(aln_set) == 1) return(result)
    
    if (length(unique(seqnames(aln_set))) > 1) {     
      new_name <- "chimera:different_chrs"
      if (keep_multichr == FALSE){
        result <- rep(FALSE, length(aln_set))
      }

    } else if (length(unique(strand(aln_set))) > 1) {
      # Case: inversion (possibly overlapping)
      # Note - haven't considered inversion + something else
      new_name <- "chimera:inversion"
    
    } else {
      # Note: alignments are sorted
      # (start ranges 2 onwards) - (end ranges 1 - second last):
      gap_lengths <- start(aln_set)[-1] - (end(aln_set)[1:(length(aln_set)-1)] - 1)

      if (all(gap_lengths < 0)){ new_name <- "chimera:duplication"
      }else if (all(gap_lengths > 0)){ new_name <- "chimera:long_gap"
      }else new_name <- "chimera:complex"
    }
    
    read_types[names(aln_set)[1]] <<- new_name
    return(result)
  },
  
  mergeLongGaps = function(){
    
    #___________________________________
    # To do:
    
    # Remove the clipped ranges from regions to be merged
    #qranges <- unlist(cigarRangesAlongQuerySpace(cigar(alns)))
    #is_aligned <- which(!unlist(explodeCigarOps(cigar(alns))) %in% c("H", "S"))
    #aln_ranges <- qranges[is_aligned]
    #rranges <- cigarRangesAlongReferenceSpace(cigar(alns))   
    
    #start_ranges <- min(which(width(rranges) > 0)) # The first aligned range
    #start_ranges[1] <- 1 # Or the first range if it's the first section 
    #end_ranges <- sapply(width(rranges), function(x) max(which(x > 0))) # The last aligned,
    ## or the last if it's the last aligned range,
    #end_ranges[length(end_ranges)] <- length(rranges[[length(rranges)]])
    
     #___________________________________
   
  },
  
  readsToTarget = function(alns, target_start, target_end, rc){
    # Narrow the reads and cigar strings to the target region, set these attributes
    
    # Narrowing example:
    # 3-4-5-6-7-8-9-10 Read
    #     5-6-7-8      Target sequence
    # target_start 5 - (read_start 3 - 1) = index 3
    # target_end 8 - target_start 5 + cigstart 3 = index 6
    # Note that alignments and cigars are reversed if strand is -ve,
    #  but start is still genomic
        
    clip_starts <- rep(0, length(alns))
    is_clipped <- which(sapply(.self$cigar_ops, '[[', 1) == "S")
    clip_starts[is_clipped] <- unlist(lapply(width(.self$query_ranges[is_clipped]), '[[', 1))
    cig_starts <- target_start - (start(alns) - 1)
    cig_ends <- target_end - target_start  + cig_starts  
    
    locs <- .self$findDeletions(cig_starts, cig_ends)
    temp <- cigarNarrow(cigar(alns), locs$starts, locs$ends)  
    new_starts <-  attr(temp, "rshift") + 1 + clip_starts 
    # + 1 as rshift is number removed not starting point
    
    # new_starts are genomic offsets wrt current genomic starts.
    # To get the corresponding alignments, translate the genomic offsets into read locs 
    nranges_cut <- max(which(start(.self$ref_ranges) <= locs$start)) 
    
    # Note that if there is an insertion at the first position of the target site,
    # cigarNarrow removes it
    for( i in seq_along(.self$cigar_ops)){
        prior_ins <- grep("I", .self$cigar_ops[[i]][1:nranges_cut[i]])
        
        if (length(prior_ins) >= 1){    
          new_starts[i] <- new_starts[i] + sum(width(.self$query_ranges[[i]][prior_ins]))
        }
        
        prior_dels <- grep("D", .self$cigar_ops[[i]][1:nranges_cut[i]])
        if (length(prior_dels) >= 1){
          new_starts[i] <- new_starts[i] - sum(width(.self$ref_ranges[[i]][prior_dels]))
        }  
    }
    
    # Set cigars and update ranges attributes
    cigs <- as.character(temp)
    query_ranges <<- cigarRangesAlongQuerySpace(cigs)
    shift_starts <- start(alns) + attr(temp, "rshift") -1
    genome_ranges <<- shift(cigarRangesAlongReferenceSpace(cigs), shift_starts)
        
    seq_lens <- sapply(.self$query_ranges, function(x) sum(width(x)))   
    seqs <- subseq(mcols(alns)$seq, start = new_starts, width = seq_lens)
    
    if (rc == TRUE){
       cigs <- unname(sapply(cigs, reverseCigar))
       query_ranges <<- cigarRangesAlongQuerySpace(cigs)
       # Note that even if strand is -ve, it is displayed wrt reference strand in bam
       seqs <- reverseComplement(seqs)
       genome_ranges <<- IRangesList(lapply(.self$genome_ranges, rev))
    }
    
    ref_ranges <<- cigarRangesAlongReferenceSpace(cigs)

    cigar_ops <<- explodeCigarOps(cigs)
    alns <<- GAlignments(seqnames = seqnames(alns), pos=start(alns) + attr(temp, "rshift"),
                 cigar = cigs, names = names(alns), strand = strand(alns),
                 seqlengths = seqlengths(alns), seq = seqs)
    
  },
  
  findDeletions = function(starts_wrt_read, ends_wrt_read){
      # Get coordinates for narrowing cigars.  For reads with a deletion
      # spanning one or both ends of the target location, narrow the 
      # cigar to encompass the deletion
      #
      # Deletions may be coded as either "D" or "N" (splice junction), 
      # depending upon the mapping software
            
      is_del <- rep(0, length(starts_wrt_read))
      del_chars <- c("N", "D")
      
      del_locs <- lapply(seq_along(starts_wrt_read), function(i){
          rs <- starts_wrt_read[i]
          re <- ends_wrt_read[i]
          fo <- findOverlaps(.self$ref_ranges[[i]], IRanges(rs, re))@queryHits
          fo
      })
      
      target_ops <- lapply(seq_along(del_locs), function(i){   
          .self$cigar_ops[[i]][del_locs[[i]]]
      })

      start_del <- sapply(target_ops, function(x) x[1] %in% del_chars)
      end_del <- sapply(target_ops, function(x) tail(x, n = 1) %in% del_chars)
      dels <- cbind(start_del, end_del)
            
      new_starts <- rep(0, length(starts_wrt_read)) 
      new_ends <- rep(0, length(ends_wrt_read))    
      for (i in 1:nrow(dels)){
          r <- unname(dels[i,])
          
          if (identical(r,c(FALSE, FALSE))){
              is_del[i] <- 0
              # return original starts
              new_starts[i] <- starts_wrt_read[i]
              new_ends[i] <- ends_wrt_read[i]
              
          }else if (identical(r,c(TRUE, TRUE))){
              # return new start points
              is_del[i] <- 1
              temp <- reduce(.self$ref_ranges[[i]][del_locs[[i]]])    
              new_starts[i] <- start(temp) -1 
              new_ends[i] <- end(temp) + 1
              
          }else if (identical(r,c(TRUE, FALSE))){
              is_del[i] <- 1
              # return new starting point, old end
              new_starts[i] <- start(.self$ref_ranges[[i]][del_locs[[i]][1]]) -1
              new_ends[i] <- ends_wrt_read[i]
              
          }else if (identical(r,c(FALSE, TRUE))){
              # return new end, old start
              is_del[i] <- 1
              new_starts[i] <- starts_wrt_read[i]
              temp <- del_locs[[i]]
              new_ends[i] <- end(.self$ref_ranges[[i]][temp[length(temp)]]) + 1
          }
      }
      long_dels <<- is_del
      return(list(starts = new_starts, ends = new_ends)) 
  },
  
  getInsertionSeqs = function(){ 
    # Note that the start of a ref_ranges insertion is its genomic end (rightmost base)
    
    ins <- lapply(.self$cigar_ops, function(x) which(x == "I"))
    idxs <- rep(1:length(ins), lapply(ins, length))
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
    dels <- lapply(.self$cigar_ops, function(x) which(x %in% c("D", "N")))       
    del_idxs <- rep(1:length(dels), lapply(dels, length))
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

    vars
  },
  
  getShortCigars = function(match_string = "No mutation", genome_to_target = NULL, rc = FALSE){
    # start_to_target, if provided, is a vector with names being genomic start coords
    # and values being coords with respect to the target site  
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
       # if (mutn["long_del"] == 1){
       #   sprintf('%s-%s:D', as.numeric(mutn["start"]),  as.numeric(mutn["end"]))
       # }else if (mutn["op"] %in% c("N", "D")){ 
       #   sprintf('%s:%sD', as.numeric(mutn["start"]), as.numeric(mutn["rwidth"]))
       # }else if (mutn["op"] == "I"){
       #   sprintf('%s:%sI', as.numeric(mutn["start"]), as.numeric(mutn["qwidth"]))
       # }
       
       #if (mutn["long_del"] == 1){
       #  sprintf('(%s,%s):D', as.numeric(mutn["start"]),  as.numeric(mutn["end"]))
       #}else if (mutn["op"] %in% c("N", "D")){ 
       #  sprintf('%s:%sD', as.numeric(mutn["start"]), as.numeric(mutn["rwidth"]))
       #}else if (mutn["op"] == "I"){
       #  sprintf('%s:%sI', as.numeric(mutn["start"]), as.numeric(mutn["qwidth"]))
       #}
       if (mutn["op"] %in% c("N", "D")){ 
         sprintf('%s:%sD', as.numeric(mutn["start"]), as.numeric(mutn["rwidth"]))
       }else if (mutn["op"] == "I"){
         sprintf('%s:%sI', as.numeric(mutn["start"]), as.numeric(mutn["qwidth"]))
       }
       
    }
    result <- rep(match_string, length(idxs))
    if (nrow(info) == 0) return(result)
    short_ops <- apply(info, 1, getShortOp)
    short_ops <- split(short_ops, cig_idxs)
    short_ops <- sapply(short_ops, function(x) do.call(paste, c(as.list(x),sep = ",")))
    result[as.numeric(names(short_ops))] <- short_ops
    result
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
             "genome_to_target")
)

CrisprSet$methods(
  initialize = function(bam_fnames, ref, target_chr, target_start, target_end, rc = FALSE,
                        short_cigars = FALSE, names = NULL, exclude_ranges = GRanges(), 
                        exclude_names = NULL, merge_chimeras = TRUE, 
                        renumbered = FALSE, target_loc = NA, ...){
    
    # TO DO - DECIDE WHETHER TO KEEP EMPTY RUNS
    target <<- GRanges(target_chr, IRanges(target_start, target_end))
    
    ref <<- ref 
    crispr_runs <<- lapply(bam_fnames, CrisprRun$new, target_chr, target_start, 
                           target_end, rc = rc, exclude_ranges = exclude_ranges,
                           exclude_names = exclude_names, merge_chimeras = merge_chimeras,
                           ...) 
    
    
    if (! is.null(names)) names(.self$crispr_runs) <- names
    nonempty_runs <<-  sapply(.self$crispr_runs, function(x) {
                              ! class(x$alns) == "uninitializedField"})
    
    .self$crispr_runs <<- .self$crispr_runs[.self$nonempty_runs]
    .self$countCigars(short_cigars, renumbered, target_loc, target_start, target_end, rc)
  },
   
  countCigars = function(short = FALSE, renumbered = FALSE, target_loc = NA, target_start = NA,
                         target_end = NA, rc = FALSE, ...){
    
    # Note that this function does not consider starts, two alignments starting at
    # different locations but sharing a cigar string are considered equal
    
    if (short == TRUE){
      g_to_t = NA
      if (renumbered == TRUE){
        if (any(is.na(c(target_loc, target_start, target_end, rc)))){
          stop("Must specify target_loc (cut site), target_start, target_end and rc
                for renumbering")
        }
        g_to_t <- genomeToTargetLocs(target_loc, target_start, target_end, rc)
      }
          
      cig.by.run <- lapply(.self$crispr_runs,#[.self$nonempty_runs], 
                    function(crun) crun$getShortCigars(genome_to_target = g_to_t, ...))
    }else {
      cig.by.run <- lapply(.self$crispr_runs,#[.self$nonempty_runs], 
                            function(crun) cigar(crun$alns))
    }
    unique_cigars <- unique(unlist(cig.by.run))
    m <- as.matrix(sapply(cig.by.run, function(x) table(x)[unique_cigars]))
    m[is.na(m)] <- 0
    rownames(m) <- unique_cigars
    m <- m[order(rowSums(m), decreasing = TRUE),]
    cigar_freqs <<- m
  },
  
  heatmapCigarFreqs = function(as_percent = FALSE, x_size = 10, x_axis_title = NULL,
                               x_angle = 90, annotate_counts = TRUE){
    
    p <- cigarFrequencyHeatmap(.self$cigar_freqs, as_percent, x_size, x_axis_title,
                               x_angle, annotate_counts)
    return(p)
  },
  
  plotVariants = function(freq_cutoff = 0, top_n = nrow(.self$cigar_freqs), ...){
    # var_freq_cutoff = i (integer) only plot variants that occur >= i times
    # top_n = total number of variants to plot
    # ... arguments for plotAlignments
    
    upto <- min(top_n, max(which(rowSums(.self$cigar_freqs) >= freq_cutoff)))
    cig_freqs <- .self$cigar_freqs[1:upto,]
    alns <- .self$makePairwiseAlns(cig_freqs)
    if (class(.self$insertion_sites) == "uninitializedField" | 
        !("cigar" %in% colnames(.self$insertion_sites))){
      .self$getInsertions() 
    }
    p <- plotAlignments(.self$ref, alns, .self$insertion_sites, ...)    
    return(p)
  },
  
  getInsertions = function(with_cigars = TRUE){
    if (with_cigars == FALSE){
      all_ins <- do.call(rbind, lapply(.self$crispr_runs, function(x) x$insertions))
    } else {
      all_ins <- do.call(rbind, lapply(.self$crispr_runs, function(x) {
                ik <- x$ins_key 
                v <- data.frame(ik, cigar(x$alns[as.integer(names(ik))]))
                v <- v[!duplicated(v),]
                v <- v[order(v$ik),]
                cbind(x$insertions[v[,1],], cigar = v[,2])
              }))
    }
    insertion_sites <<- all_ins[order(all_ins$start, all_ins$seq),]
  },
  
  makePairwiseAlns = function(cig_freqs = .self$cigar_freqs, ...){
    # Get alignments by cigar string, make the alignment for the consensus

    cigs <- unlist(lapply(.self$crispr_runs, function(x) cigar(x$alns)), use.names = FALSE)
    splits <- split(seq_along(cigs), cigs)  
    splits <- splits[match(rownames(cig_freqs), names(splits))]
    
    x <- lapply(.self$crispr_runs, function(x) x$alns)
    all_alns <- do.call(c, unlist(x, use.names = FALSE))
  
    seqs <- c()
    starts <- c()
  
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
   alns
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
    
    # TO DO - MUST BE A NICER WAY OF DOING THIS!
    gs <- min(unlist(lapply(.self$crispr_runs, function(x) start(unlist(x$genome_ranges)))))
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


