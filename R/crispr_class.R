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
  # THIS FUNCTION NOT WORKING BECAUSE OF VIEW MEANS
  # Chimeras: GAlignments obj
  chimera_cov <- coverage(chimeras)
  sl <- slice(chimera_cov, lower = min_cov)
  slgr <- as(sl, "GRanges")
  vm <- viewMeans(sl)
  mcols(slgr)$avg <- unlist(vm)
  return(slgr)
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


getCodonFrame <- function(txdb, target_chr, target_start, target_end){
  # Target_chr must match txdb
  require(VariantAnnotation)
  refLocsToLocalLocs(GRanges(target_chr, IRanges(target_start, target_end)), txdb)
  
  # check that all transcripts have the same frame
  
}

nucleotideToAA <- function(seqs, txdb, target_start, target_end){
    
}

