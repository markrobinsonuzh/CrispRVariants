#'@title Get variant counts
#'@description Returns a matrix of counts where rows are sequence variants an columns are samples       
#'@param obj An object containing variant counts
#'@param ... Additional arguments
#'@return A matrix of counts where rows are variants and columns are samples
#'@author Helen Lindsay
#'@rdname variantCounts
#'@export
setGeneric("variantCounts", function(obj, ...) {
  standardGeneric("variantCounts")})

#'@param freq.cutoff (Integer n) Return variants with total count greater than or 
#'equal to n
#'@param top.n  (Integer n) If specified, return variants ranked at least n according
#' to frequency across all samples (Default: 0, i.e. no cutoff)
#'@param include.chimeras Should chimeric reads be included in the counts table?
#'(Default: TRUE)
#'@rdname variantCounts
#'@examples
#'data("gol_clutch1")
#'
#'#Return a matrix of the 5 most frequent variants
#'variantCounts(gol, top.n = 5)
setMethod("variantCounts", signature("CrisprSet"),
          function(obj, ..., top.n = NULL, freq.cutoff = 1, 
                   include.chimeras = TRUE){
    
    if (is.null(top.n) & freq.cutoff == 0){
        return(obj$.getFilteredCigarTable(include.chimeras = include.chimeras))
    }
    
    top.n <- ifelse(is.null(top.n), nrow(obj$cigar_freqs), top.n)
    return(obj$.getFilteredCigarTable(top.n, freq.cutoff, include.chimeras))
})
          

#'@title Get mutation efficiency
#'@description  Returns the percentage of sequences that contain at least one mutation.       
#'@param obj An object containing variant counts
#'@param ... additional arguments
#'@author Helen Lindsay
#'@rdname mutationEfficiency
#'@export
setGeneric("mutationEfficiency", function(obj, ...) {
  standardGeneric("mutationEfficiency")})
               
#'@param snv  Single nucleotide variants (SNVs) may be considered as mutations ("include"),
#'treated as ambiguous sequences and not counted at all ("exclude"), or treated as 
#'non-mutations, e.g. sequencing errors or pre-existing SNVs ("non_variant")
#'@param include.chimeras Should chimeric alignments be counted as variants 
#'when calculating mutation efficiency (Default: TRUE
#'@param exclude.cols A vector of names or indices of columns in the variant counts table
#'that will not be considered when counting mutation efficiency  
#'@rdname mutationEfficiency
#'@return A vector of efficiency statistics per sample and overall
#'@examples
#'data("gol_clutch1")
#'mutationEfficiency(gol)
setMethod("mutationEfficiency", signature("CrisprSet"),
          function(obj, ..., snv = c("include","exclude","non_variant"),
                   include.chimeras = TRUE, exclude.cols = NULL){
    return(obj$mutationEfficiency(snv = snv, chimeras = chimeras, 
                                  exclude_cols = exclude.cols))
})
