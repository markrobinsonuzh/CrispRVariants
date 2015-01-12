
# When given nothing, returns the variant frequency table
# When given labels, 

# groupVariants

#'@title Counts variants within a target region 
setGeneric("variantCounts", function(...) {
  standardGeneric("variantCounts")})

#setMethod("variantCounts", signature("CrisprSet"),
#          function(obj, ... ){})
            
# heatmapCigarFreqs
# plotVariants
# plotFrequencySpectrum
# mutationEfficiency
# getVariantTable


#'@title Plot alignments with respect to a reference sequence  
#'@rdname plotAlignments
#'@export
setGeneric("plotAlignments", function(obj, ...) {
  standardGeneric("plotAlignments")})

#'@rdname plotAlignments
#'@description Wrapper for CrisprSet$plotVariants.  Optionally filters a 
#'CrisprSet frequency table, then plots variants.  More information in
#'\code{\link[crispRvariants]{CrisprSet}} 
#'@param freq_cutoff i (integer) only plot variants that occur >= i times
#' (default: 0, i.e no frequency cutoff)
#'@param top_n (integer) Plot only the n most frequent variants 
#' (default: plot all)
#'@param renumbered If TRUE, the x-axis is numbered with respect to the target
#' (default: TRUE)
setMethod("plotAlignments", signature("CrisprSet"),  
          function(obj, ..., freq_cutoff = 0, 
                   top_n = nrow(obj$cigar_freqs),
                   renumbered = obj$pars["renumbered"]) {
  plot_obj <- obj$plotVariants(freq_cutoff = freq_cutoff, top_n = top_n, 
                               renumbered = renumbered, ...)
  
  return(plot_obj)
})
