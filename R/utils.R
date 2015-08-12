#'@title Merge two CrisprSets
#'@description Merge two CrisprSet objects sharing a reference and target location
#'@author Helen Lindsay
#'@rdname mergeCrisprSets
setGeneric("mergeCrisprSets", function(x,y, ...) {
  standardGeneric("mergeCrisprSets")})


#'@rdname mergeCrisprSets
#'@param x A CrisprSet object
#'@param y A second CrisprSet object
#'@param order A list of sample names, matching the names in x and y,
#'specifying the order of the samples in the new CrisprSet. (Not implemented yet)
#'@param ... extra arguments
#'@return A merged CrisprSet object
#'@examples
#'# Load the metadata table
#'md_fname <- system.file("extdata", "gol_F1_metadata_small.txt", package = "crispRvariants")
#'md <- read.table(md_fname, sep = "\t", stringsAsFactors = FALSE)
#'
#'# Get bam filenames and their full paths
#'bam_fnames <- sapply(md$bam.filename, function(fn){
#'  system.file("extdata", fn, package = "crispRvariants")})
#'
#'reference <- DNAString("GGTCTCTCGCAGGATGTTGCTGG")
#'gd <- GRanges("18", IRanges(4647377, 4647399), strand = "+")
#'
#'crispr_set1 <- readsToTarget(bam_fnames[c(1:4)], target = gd, 
#'       reference = reference, names = md$experiment.name[1:4], target.loc = 17)
#'crispr_set2 <- readsToTarget(bam_fnames[c(5:8)], target = gd, 
#'       reference = reference, names = md$experiment.name[5:8], target.loc = 17)
#'mergeCrisprSets(crispr_set1,crispr_set2)
#'@export
setMethod("mergeCrisprSets", signature(x = "CrisprSet", y = "CrisprSet"),
          function(x,y, ..., order = NULL){
            # To do: add order vector
            
            warn_msg <- "CrisprSets must have the same %s for merging"
            ref = x$ref
            if (ref != y$ref){
              stop(sprintf(warn_msg, "reference sequence"))
            }
            target = x$target
            if (target != y$target){
              stop(sprintf(warn_msg, "guide sequence"))
            }
            t.loc = x$pars$target.loc
            if (t.loc != y$pars$target.loc){
              stop(sprintf(warn_msg, "target location (zero point)"))
            }                       
            names <- c(names(x), names(y))
            
            cruns <- c(x$crispr_runs, y$crispr_runs)
            cset <- CrisprSet(cruns, reference = ref, target = target, 
                      names = names, target.loc = t.loc)
            return(cset)
          })




#'@title Count the number of reads containing an insertion or deletion
#'@description Counts the number of reads containing a deletion or insertion 
#'(indel) of any size in a set of aligned reads.  
#'For countDeletions and countInsertions Reads may be filtered according to 
#'whether they contain more than one indel of the same or different types.
#'@author Helen Lindsay
#'@param alns The aligned reads
#'@param multi.del  If TRUE, returns the exact number of deletions,
#'i.e., if one read contains 2 deletions, it contributes 2 to the
#'total count (default: FALSE)
#'@param del.and.ins If TRUE, counts deletions regardless of whether 
#'reads also contain insertions.  If FALSE, counts reads that contain 
#'deletions but not insertions (default: FALSE)
#'@param del.ops Cigar operations counted as deletions.  Default: c("D")
#'@param ... extra arguments
#'@rdname indelCounts
#'@examples
#'bam_fname <- system.file("extdata", "gol_F1_clutch_2_embryo_4_s.bam",
#'                          package = "crispRvariants")
#'bam <- GenomicAlignments::readGAlignments(bam_fname, use.names = TRUE)
#'countDeletions(bam)
#'countInsertions(bam)
#'countIndels(bam)
#'indelPercent(bam)
#'@export
setGeneric("countDeletions", function(alns, ...) {
  standardGeneric("countDeletions")})


#'@return countDeletions: The number of reads containing a deletion (integer)
#'@rdname indelCounts
setMethod("countDeletions", signature("GAlignments"),
          function(alns, ..., multi.del = FALSE, del.and.ins = FALSE,
                   del.ops=c("D")){
 
  cigar_ops <- CharacterList(explodeCigarOps(cigar(alns)))
  
  if (multi.del == TRUE){
    has_del <- any(cigar_ops %in% del.ops)
  } else{
    has_single_del <- sum(cigar_ops %in% del.ops) == 1
  }
  
  if (del.and.ins){
    if (multi.del)  return(sum(has_del))
    return(sum(has_single_del))      
  }
  
  has_ins <- any(cigar_ops == "I")
  
  if (multi.del){
    # multi.del == TRUE and del.and.ins == FALSE
    return(sum(has_del & ! has_ins))
  }
  
  # multi.del == FALSE and del.and.ins == FALSE
  return(sum(has_single_del & ! has_ins))
})



#'@rdname indelCounts
#'@export
setGeneric("countInsertions", function(alns, ...) {
  standardGeneric("countInsertions")})


#'@param multi.ins  If TRUE, returns the exact number of insertions,
#'i.e., if one read contains 2 insertions, it contributes 2 to the
#'total count (default: FALSE)
#'@param ins.and.del If TRUE, counts insertions regardless of whether 
#'reads also contain deletions  If FALSE, counts reads that contain 
#'insertions but not deletions (default: FALSE) 
#'@return countInsertions: The number of reads containing an insertion (integer)
#'@rdname indelCounts
setMethod("countInsertions", signature("GAlignments"),
          function(alns, ..., ins.and.del = FALSE, multi.ins = FALSE, del.ops = c("D")){
  
  cigar_ops <- CharacterList(explodeCigarOps(cigar(alns)))
  
  if (multi.ins == TRUE){
    has_ins <- any(cigar_ops == "I")
  } else{
    has_single_ins <- sum(cigar_ops == "I") == 1
  }
  
  if (ins.and.del){
    if (multi.ins)  return(sum(has_ins))
    return(sum(has_single_ins))
  }
  
  has_del <- any(cigar_ops %in% del.ops)
  
  if (multi.ins) return(sum(has_ins &! has_del))
  return(sum(has_single_ins &! has_del))
})


#'@rdname indelCounts
#'@export
setGeneric("countIndels", function(alns) {
  standardGeneric("countIndels")})

#'@return countIndels: The number of reads containing at least one insertion 
#'@rdname indelCounts
setMethod("countIndels", signature("GAlignments"),
          function(alns){
  cigar_ops <- CharacterList(explodeCigarOps(cigar(alns)))
  return(sum(any(cigar_ops %in% c("I", "D", "N"))))
})

#'@rdname indelCounts
setGeneric("indelPercent", function(alns) {
  standardGeneric("indelPercent")})


#'@rdname indelCounts
#'@return indelPercent: The percentage of reads containing an insertion or 
#'deletion (numeric)
#'@export
setMethod("indelPercent", signature("GAlignments"),
          function(alns){
  return((countIndels(alns) / length(alns))*100)
})
