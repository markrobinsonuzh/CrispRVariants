library(GenomicAlignments)
library(Rsamtools)
library(parallel)
library(ggplot2)
library(reshape2)
#library(sangerseqR)

# To do - deal with soft clipping in the target region
# To do - check are empty alignments dealt with correctly?
# To do - renumber wrt target
# To do - check for sapply with potential for unwanted simiplification
# To do - read bam - exclude should be granges?



writeFastq <- function(outf, vals){
    o <- file(outf, "a")
    seqname <- sprintf("@%s", vals$seqname)
    qualname <- sprintf("+%s", vals$seqname)
    writeLines(c(seqname, vals$seq, qualname, vals$quals) , o)
    close(o)
}

abifToTrimmedFastq <- function(seqname, fname, outfname, trim = TRUE, cutoff = 0.05, 
                               min_seq_len = 20){
    # Translation of python function
    abif <- read.abif(fname)
    if (is.null(abif@data$PCON.2)){
       print(sprintf("failed on %s", seqname))
       return()
    }
    
    # Hack to remove the extra character if it exists
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
                "quals" = rawToChar(as.raw(num_quals[trim_start:trim_finish]+33))))
    return()
}

mergeChimericAlns <- function(alns){
  # Chimeras must have multiple seqs on the same chromosome and strand
  if (length(alns) == 1) return(alns)
  if (length(unique(seqnames(alns))) > 1) {
    print(c("Ignored chimera from different chromosomes", seqnames(alns)))
    return(alns)
  }
  if (length(unique(strand(alns))) > 1) {
    print(c("Ignored inversion", seqnames(alns))
    return(alns)
  }
     
  qranges <- unlist(cigarRangesAlongQuerySpace(cigar(alns)))
  is_aligned <- which(!unlist(explodeCigarOps(cigar(alns))) %in% c("H", "S"))
  aln_ranges <- qranges[is_aligned]
  
  # Aligned regions must not overlap
  if (! length(findOverlaps(qranges[ops], ignoreSelf = TRUE)) == 0) return(alns)
    
  # Merge the cigars.  
  # Is it safe to assume that the original sequence is not hard clipped?    
  rranges <- cigarRangesAlongReferenceSpace(cigar(alns))   
  
  start_ranges <- min(which(width(rranges) > 0)) # The first aligned range
  start_ranges[1] <- 1 # Or the first range if it's the first section 
  end_ranges <- which.max(width(rranges)) # The last aligned range, or the last if it's
  end_ranges[length(end_ranges)] <- length(rranges[[length(rranges)]]) # the last section 
  gap_lengths <- end(alns)[1:length(alns)-1] - (start(alns)[-1] - 1)
  
    

}


readBam <- function(bam_fname, target_chr, target_start, target_end, merge_chimeras = FALSE,
                    exclude = GRanges()){
  # This function reads a bam file and classifies reads as on_target if they 
  # completely span the target, region, and off_target otherwise. 
  #
  # Note that unmapped reads are discarded by readGAlignments, 
  # and are not counted as "off-target" reads
  #
  # "exclude" should be a GRanges object of regions that should not be counted,
  #  e.g. primer or cloning vector sequences that have a match in the genome   
  #
  # If merge_chimeras = TRUE, chimeric reads (alignments with a large gap, split 
  # over two lines) are merged if both sequences are from the same chromosome, 
  # do not overlap, and are aligned to the same strand.
  # It is assumed that sequences with two alignments are chimeras, not alternate mappings


  param <- ScanBamParam(what = c("seq"))
  bam <- readGAlignments(bam_fname, param = param, use.names = TRUE)
  if (length(bam) == 0) return(list("on_target" = NULL, "off_target" = NULL))
  
  bam <- bam[setdiff(1:length(bam), findOverlaps(bam, exclude)@queryHits)]
      
  if (merge_chimeras == TRUE){
    # Split the bam file by seq name, preserving the original order
    bam_by_name <-split(bam, factor(names(bam), levels = unique(names(bam))))
    result <- sapply(bam_by_name, mergeChimericAlns)
  }  
  is_on_target <- factor(as.character(start(bam) <= target_start & end(bam) >= target_end & 
                   seqnames(bam) == target_chr), levels = c("TRUE", "FALSE"))
  result <- split(bam, is_on_target)
  names(result) <- c("on_target", "off_target")
  result
}   

reverseCigar <- function(cigar){
   cigar.widths <- rev(strsplit(cigar, '[A-Z]')[[1]])
   cigar.ops <- rev(explodeCigarOps(cigar)[[1]])
   paste0(cigar.widths,cigar.ops, collapse = "")
}

seqsToAln <- function(cigar, dnaseq, del_char = "-"){
    wrt_ref <- cigarRangesAlongReferenceSpace(cigar)[[1]]
    wrt_read <- cigarRangesAlongQuerySpace(cigar)[[1]]
    ops <- explodeCigarOps(cigar)[[1]]
    segs <- as.character(Views(dnaseq, wrt_read))
    segs[which(ops == "I")] <- ""
    for (j in which(ops == "D")){
         segs[j] <- paste0(rep(del_char, width(wrt_ref[j])), collapse = "")
      }
    result <- paste0(segs, collapse = "")
    result
}

seqToGappedAln <- function(cigar, dnaseq, gap_locs){

}


#_______________________________________________________________________________________

CrisprRun = setRefClass(
  Class = "CrisprRun",
  fields = c("alns",
             "name",
             "query_ranges",
             "ref_ranges",
             "cigar_ops",
             "off_target",
             "long_dels",
             "insertions")
)

CrisprRun$methods(
  initialize = function(bam_fname, target_chr, target_start, target_end, rc = FALSE, 
                        name = NULL){
    
      name <<- ifelse(is.null(name), bam_fname, name)
      
      galns <- readBam(bam_fname, target_chr, target_start, target_end)
      off_target <<- galns[["off_target"]]
      on_target <- galns[["on_target"]]
      if (is.null(on_target)) return()

      ref_ranges <<- cigarRangesAlongReferenceSpace(cigar(on_target))
      query_ranges <<- cigarRangesAlongQuerySpace(cigar(on_target))
      cigar_ops <<- explodeCigarOps(cigar(on_target)) 
      .self$readsToTarget(on_target, target_start, target_end, rc)
      .self$getInsertionSeqs()
  },
  
  readsToTarget = function(alns, target_start, target_end, rc){
    # Narrow the reads and cigar strings to the target region, set these attributes
    
    # Narrowing example:
    # 3-4-5-6-7-8-9-10 Read
    #     5-6-7-8      Target sequence
    # target_start 5 - (read_start 3 - 1) = index 3
    # target_end 8 - target_start 5 + cigstart 3 = index 6
    
    clip_starts <- rep(0, length(alns))
    is_clipped <- which(sapply(.self$cigar_ops, '[[', 1) == "S")
    clip_starts[is_clipped] <- sapply(width(.self$query_ranges[is_clipped]), '[[', 1)
    
    cig_starts <- target_start - (start(alns) - 1)
    cig_ends <- target_end - target_start  + cig_starts
    locs <- .self$findDeletions(cig_starts, cig_ends)

    temp <- cigarNarrow(cigar(alns), locs$starts, locs$ends)    
    new_starts <-  attr(temp, "rshift") + 1 + clip_starts 
    # + 1 as rshift is number removed not starting point
    
    # Check for insertions prior to the new start sites, these affect the
    # new starting locations with respect to the read sequences
    cut_ranges <- sapply(.self$query_ranges[start(.self$query_ranges) < new_starts], length)
    for( i in seq_along(.self$cigar_ops)){
        prior_ins <- grep("I", .self$cigar_ops[[i]][1:cut_ranges[i]])
        if (length(prior_ins) >= 1){    
          new_starts[i] <- new_starts[i] + sum(width(.self$query_ranges[[i]][prior_ins]))
        }
    }
    
    # Set cigars and update ranges attributes
    cigs <- as.character(temp)
    query_ranges <<- cigarRangesAlongQuerySpace(cigs)
    seq_lens <- sapply(.self$query_ranges, function(x) sum(width(x)))   
    seqs <- subseq(mcols(alns)$seq, start = new_starts, width = seq_lens)
      
    if (rc == TRUE){
       cigs <- unname(sapply(cigs, reverseCigar))
       query_ranges <<- cigarRangesAlongQuerySpace(cigs)
       seqs <- reverseComplement(seqs)
    }
    
    ref_ranges <<- cigarRangesAlongReferenceSpace(cigs)
    cigar_ops <<- explodeCigarOps(cigs)

    alns <<- GAlignments(seqnames = seqnames(alns), pos=attr(temp, "rshift"),
                 cigar = cigs, names = names(alns), strand = strand(alns),
                 seqlengths = seqlengths(alns), seq = seqs)
    
  },
  
  findDeletions = function(starts_wrt_read, ends_wrt_read){
      # Deletions may be coded as either "D" or "N" (splice junction), 
      # depending upon the mapper
      
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
    ins <- sapply(.self$cigar_ops, function(x) which(x == "I"))
    tseqs <- as.character(mcols(.self$alns)$seq)[sapply(ins, length) > 0]
    if (length(tseqs) == 0) {
        insertions <<- data.frame(start = character(), seq = character(), count = character())
        return()
    }    
    qranges <- unlist(.self$query_ranges[ins])
    ins_seqs <- sapply(seq_along(tseqs), function(i) as.character(Views(tseqs[i],qranges[i])))    
    ins_starts <- start(unlist(.self$ref_ranges[ins]))
    df <- data.frame(start = ins_starts, seq = ins_seqs)
    df$seq <- as.character(df$seq)
    insertions <<- aggregate(rep(1, nrow(df)), by = as.list(df), FUN = table)
    colnames(.self$insertions) <- c("start", "seq", "count")
  },
  
  renumberCigars = function(cut_base, relative_locs, genomic_locs = NULL, genomic = FALSE){
      # cut_base may be either relative to the start of the target sequence, or genomic
      # cut_base is the left side of the cut site
      # (Will be numbered -1)
      
      # Example:  cut_base = 5
      # Before: 1  2  3  4  5  6  7  8 
      # After: -5 -4 -3 -2 -1  1  2  3
      # Left =  original - cut_base - 1
      # Right = original - cut_base
      
      if (genomic == TRUE){
        result <- genomic_locs
        left <- which(genomic_locs <= cut_base)
        result[locs <= cut_base] 
      }
      
  },
  
  getShortCigars = function(match_string = "No mutation"){
      idxs <- sapply(.self$cigar_ops, function(x) which(x != "M"))
      ops <- sapply(seq_along(idxs), function(i) .self$cigar_ops[[i]][idxs[[i]]])
      rranges <- .self$ref_ranges[idxs]
      qranges <- .self$query_ranges[idxs]
  
      cig_idxs <- unlist(sapply(seq_along(ops), function(i) rep(i, length(ops[[i]]))))      
      
      info <- data.frame(op = unlist(ops), qwidth = unlist(width(qranges)), 
                start = unlist(start(rranges)), end = unlist(end(rranges)), 
                rwidth = unlist(width(rranges)), long_del = .self$long_dels[cig_idxs])
  
      getShortOp <- function(mutn){
          if (mutn["long_del"] == 1){
            sprintf('%s-%s:D', as.numeric(mutn["start"]),  as.numeric(mutn["end"]))
          }else if (mutn["op"] %in% c("N", "D")){ 
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
             "cigar_freqs")
)

CrisprSet$methods(
  initialize = function(bam_fnames, ref, target_chr, target_start, target_end, rc = FALSE,
                        short_cigars = FALSE, names = NULL){
    ref <<- ref
    crispr_runs <<- sapply(bam_fnames, CrisprRun$new, target_chr, target_start, 
                           target_end, rc = rc) 
    if (! is.null(names)) names(.self$crispr_runs) <- names
    nonempty_runs <<-  sapply(.self$crispr_runs, function(x) {
                              ! class(x$alns) == "uninitializedField"})
    .self$countCigars(short_cigars)
  },
   
  countCigars = function(short = FALSE, ...){
    if (short == TRUE){
      cig.by.run <- sapply(.self$crispr_runs[.self$nonempty_runs], 
                            function(crun) crun$getShortCigars(...) )
    }else {
      cig.by.run <- sapply(.self$crispr_runs[.self$nonempty_runs], 
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
    
    counts <- melt(.self$cigar_freqs)
    colnames(counts) <- c("Cigar", "Sample","Count")
    counts$Cigar <- factor(counts$Cigar, levels = rev(levels(counts$Cigar)))
    
    if (as_percent == TRUE){
      m <- apply(.self$cigar_freqs, 2, function(x) x/sum(x))  
      m <- melt(m)
      colnames(m) <- c("Cigar", "Sample","Percentage")  
      m$Cigar <- factor(m$Cigar, levels = rev(levels(m$Cigar)))
      g <- ggplot(m, aes(x = Sample, y = Cigar, fill = Percentage)) + geom_tile()
    }
    else{
      g <- ggplot(counts, aes(x = Sample, y = Cigar, fill = Count)) + geom_tile()
    }
    if (annotate_counts == TRUE){
        g <- g + geom_text(data = counts, aes(label = Count, fill = NULL))
    }
    g <- g + ylab(NULL) + xlab(x_axis_title) + theme_bw() +
         theme(axis.text.x = element_text(size = x_size, angle = x_angle))
    g
  },
  
  plotVariants = function(){
    # TO DO - HOW TO PASS LOTS OF OPTIONAL ARGS?
    
    alns <- .self$makePairwiseAlns()
    
  },
  
  findAllInsertions = function(){
    all_ins <- lapply(.self$crispr_runs, function(x) x$insertions)
    ins_nms <- unlist(sapply(seq_along(all_ins), function(i) rep(names(all_ins[i]), nrow(all_ins[[i]]))))
    all_ins <- do.call(rbind, all_ins)
    all_ins <- cbind(all_ins, ins_nms)[not_del]
    insertion_sites <<- all_ins[order(all_ins$start, all_ins$seq),]
  },
  
  makePairwiseAlns = function(...){
    # Get alignments by cigar string, make the alignment for the consensus
    # TO DO - REPRESENTATION OF LONG DELETIONS
    
    x <- sapply(.self$crispr_runs, function(x) mcols(x$alns)$seq)
    all_seqs <- as.character(do.call(c, unlist(x, use.names = FALSE)))
    all_cigs <- unlist(sapply(.self$crispr_runs, function(x) cigar(x$alns)), use.names =FALSE)
    not_del <- unlist(lapply(.self$crispr_runs, function(x) x$long_dels), use.names = FALSE) != 1
    seq_by_cig <- split(all_seqs[not_del], all_cigs[not_del])
    
    # Order to match cigar frequencies
    seq_by_cig <- seq_by_cig[na.omit(match(rownames(.self$cigar_freqs), names(seq_by_cig)))]
    
    # TO DO - proper consensus seq, now just taking most frequent
    seqs <- sapply(seq_by_cig, function(x) {
                   tt <- table(x)
                    names(tt)[which.max(tt)]})
    alns <- mapply(seqsToAln, names(seqs), seqs, ...)
    alns 
    
  },
  
  renumberTarget = function(target_seq, target_loc){
    # Target location is between two nucleotides
    # Should reference be a list or a character vector?
    
    new_numbering <- c(seq(-1*target_loc,-1), c(1:(length(target_seq[[1]]) - target_loc)))
  }
)


