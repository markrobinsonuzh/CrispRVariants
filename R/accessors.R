
# When given nothing, returns the variant frequency table
# When given labels, 

# groupVariants


#'@title Get variant counts
#'@description Returns a matrix of counts where rows are sequence variants an columns are samples       
#'@param obj An object containing variant counts
#'@author Helen Lindsay
#'@rdname variantCounts
#'@export
setGeneric("variantCounts", function(obj, ...) {
  standardGeneric("variantCounts")})

#'@param freq.cutoff (Integer n) Return variants with total count greater than or 
#'equal to n
#'@param top.n  (Integer n) If specified, return variants ranked at least n according
#' to frequency across all samples (Default: 0, i.e. no cutoff)
#'@rdname variantCounts
setMethod("variantCounts", signature("CrisprSet"),
          function(obj, ..., top.n = NULL, freq.cutoff = 0){
    if (is.null(top.n) & freq.cutoff == 0) return(obj$cigar_freqs)
    top.n <- ifelse(is.null(top.n), nrow(obj$cigar_freqs), top.n)
    return(obj.getFilteredCigarTable(top.n, freq.cutoff))
})
          

#'@title Get mutation efficiency
#'@description  Returns the percentage of sequences that contain at least one mutation.       
#'@param obj An object containing variant counts
#'@author Helen Lindsay
#'@rdname mutationEfficiency
#'@export
setGeneric("mutationEfficiency", function(obj, ...) {
  standardGeneric("mutationEfficiency")})
               
#'@param snv  Single nucleotide variants (SNVs) may be considered as mutations ("include"),
#'treated as ambiguous sequences and not counted at all ("exclude"), or treated as 
#'non-mutations, e.g. sequencing errors or pre-existing SNVs ("non_variant")
#'@param exclude.cols A vector of names or indices of columns in the variant counts table
#'that will not be considered when counting mutation efficiency  
#'@rdname mutationEfficiency
#'@return A vector of efficiency statistics per sample and overall
setMethod("mutationEfficiency", signature("CrisprSet"),
          function(obj, ..., snv = c("include","exclude","non_variant"),
                                                 exclude.cols = NULL){
    return(obj$mutationEfficiency(snv = snv, exclude_cols = exclude.cols))
})






            
# heatmapCigarFreqs
# plotFrequencySpectrum
# getVariantTable
