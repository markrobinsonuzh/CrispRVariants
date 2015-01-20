
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
