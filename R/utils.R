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
#'specifying the order of the samples in the new CrisprSet. 
#'@param ... extra arguments
#'@export
setMethod("mergeCrisprSets", signature(x = "CrisprSet", y = "CrisprSet"),
          function(x,y, ..., order = NULL){
            # docType methods
            ## aliases merge,CrisprSet,CrisprSet-method  
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
            
            # Reorder if required.  Can be used to subset
            #if (length(order) > 1){
            #  if (! length(order) == length(names)){
            #    stop("")
            #  }
            #}
            cruns <- c(x$crispr_runs, y$crispr_runs)
            cset <- CrisprSet(cruns, reference = ref, target = target, 
                      names = names, target.loc = t.loc)
            return(cset)
          })




#'@title Count the number of reads containing a deletion
#'@description Counts the number of reads containing a deletion of any size
#'in a set of aligned reads
#'@author Helen Lindsay
#'@param alns The aligned reads
#'@param multi.del  If TRUE, returns the exact number of deletions,
#'i.e., if one read contains 2 deletions, it contributes 2 to the
#'total count (default: FALSE)
#'@param del.and.ins If TRUE, counts deletions regardless of whether 
#'reads also contain insertions.  If FALSE, counts reads that contain 
#'deletions but not insertions (default: FALSE)
#'@param del.ops Cigar operations counted as deletions.  Default: c("D")
#'@rdname countDeletions
setGeneric("countDeletions", function(alns, ...) {
  standardGeneric("countDeletions")})


#'@rdname countDeletions
setMethod("countDeletions", signature("GAlignments"),
          function(alns, multi.del = FALSE, del.and.ins = FALSE,
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


#'@title Counts the number of reads containing an insertion
#'@description Counts the number of reads containing an insertion of any size
#'in a set of aligned reads
#'@author Helen Lindsay
#'@param alns The aligned reads
#'@param multi.ins  If TRUE, returns the exact number of insertions,
#'i.e., if one read contains 2 insertions, it contributes 2 to the
#'total count (default: FALSE)
#'@param ins.and.del If TRUE, counts insertions regardless of whether 
#'reads also contain deletions  If FALSE, counts reads that contain 
#'insertions but not deletions (default: FALSE) 
#'@param del.ops Cigar operations counted as deletions.  Default: c("D")
#'@rdname countInsertions
setGeneric("countInsertions", function(alns, ...) {
  standardGeneric("countInsertions")})


#'@rdname countInsertions
setMethod("countInsertions", signature("GAlignments"),
          function(alns, ins.and.del = FALSE, multi.ins = FALSE, del.ops = c("D")){
  
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



#'@title Counts the number of reads containing an insertion and/or deletion
#'@description Counts the number of target reads that include at least one 
#'insertion or deletion, by counting cigar operations "I" (insertion),
#'"D" (deletion) and "N" (junction operation used by some aligners)
#'@author Helen Lindsay
#'@param alns The alignments 
#'@rdname countIndels
countIndels <- function(alns){
  cigar_ops <- CharacterList(explodeCigarOps(cigar(alns)))
  return(sum(any(cigar_ops %in% c("I", "D", "N"))))
}


#'@title Counts the percentage of reads containing an insertion and/or deletion
#'@description Prints the percentage of target reads that include at least one 
#'insertion or deletion
#'@author Helen Lindsay
#'@param alns The alignments 
#'@rdname indelPercent
indelPercent <- function(alns){
  return((countIndels(alns) / length(alns))*100)
}
