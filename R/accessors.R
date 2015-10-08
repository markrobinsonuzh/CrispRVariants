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

#'@param min.freq (Float n%) Return variants with frequency at least n% in at
#'least one sample (Default: 0)
#'@param min.count (Integer n) Return variants with count greater than n
#'in at least one sample (Default: 0)
#'@param top.n  (Integer n) If specified, return variants ranked at least n according
#' to frequency across all samples (Default: 0, i.e. no cutoff)
#'@param include.chimeras Should chimeric reads be included in the counts table?
#'(Default: TRUE)
#'@param include.nonindel Should sequences without indels be returned?
#'(Default:TRUE)
#'@param result Return variants as either counts ("counts", default) or
#'proportions ("proportions")
#'@rdname variantCounts
#'@examples
#'data("gol_clutch1")
#'
#'#Return a matrix of the 5 most frequent variants
#'variantCounts(gol, top.n = 5)
setMethod("variantCounts", signature("CrisprSet"),
          function(obj, ..., top.n = NULL, min.freq = 0, min.count = 1,
              include.chimeras = TRUE, include.nonindel=TRUE,
              result = "counts"){

    if (is.null(top.n) & min.freq == 0){
        return(obj$.getFilteredCigarTable(include.chimeras = include.chimeras,
                                          include.nonindel = include.nonindel))
    }

    top.n <- ifelse(is.null(top.n), nrow(obj$cigar_freqs), top.n)
    result <- obj$.getFilteredCigarTable(top.n, min.freq, min.count,
                                         include.chimeras,
                                         include.nonindel, result)
    result
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
#'non-mutations, e.g. sequencing errors or pre-existing SNVs ("non_variant", default)
#'@param include.chimeras Should chimeric alignments be counted as variants
#'when calculating mutation efficiency (Default: TRUE
#'@param exclude.cols A vector of names of columns in the variant counts table
#'that will not be considered when counting mutation efficiency
#'@param filter.vars Variants to remove before calculating efficiency.  May be either
#'a variant size, e.g. "1D", or a particular variant/variant combination, e.g. -5:3D
#'@param filter.cols A vector of control sample names.  Any variants present in the control
#'samples will be counted as non-variant, unless they also contain another indel.  Note that
#'this is not compatible with counting snvs as variants.
#'@rdname mutationEfficiency
#'@return A vector of efficiency statistics per sample and overall
#'@examples
#'data("gol_clutch1")
#'mutationEfficiency(gol)
setMethod("mutationEfficiency", signature("CrisprSet"),
          function(obj, ..., snv = c("non_variant", "include","exclude"),
                   include.chimeras = TRUE, exclude.cols = NULL,
                   filter.vars = NULL, filter.cols = NULL){
    result <- obj$mutationEfficiency(snv = snv,
                                  include.chimeras = include.chimeras,
                                  exclude.cols = exclude.cols,
                                  filter.vars = filter.vars,
                                  filter.cols = filter.cols)
    result
})


#'@title Find frequent SNVs
#'@description Find single nucleotide variants (SNVs) above a specified frequency
#'in a table of variants.
#'@param obj An object containing variant counts
#'@param ... additional arguments
#'@author Helen Lindsay
#'@rdname findSNVs
#'@export
setGeneric("findSNVs", function(obj, ...) {
  standardGeneric("findSNVs")})


#'@rdname findSNVs
#'@param freq minimum frequency snv to return (Default: 0.25)
#'@param include.chimeras include chimeric reads when calculating SNV frequencies
#'(Default: TRUE)
#'@return A vector of SNVs and their frequencies
setMethod("findSNVs", signature("CrisprSet"),
          function(obj, ..., freq = 0.25, include.chimeras = TRUE){
            result <- obj$.getSNVs(min.freq = freq,
                                  include.chimeras = include.chimeras)
            result
          })


#'@title Get chimeric alignments
#'@description Return chimeric alignments from a collection of aligned sequences
#'@param obj An object containing aligned sequences
#'@param ... additional arguments
#'@author Helen Lindsay
#'@rdname getChimeras
#'@export
setGeneric("getChimeras", function(obj, ...) {
  standardGeneric("getChimeras")})

#'@rdname getChimeras
#'@param sample The sample name or sample index to return
#'@return A GAlignment object containing the chimeric read groups
#'@examples
#'data("gol_clutch1")
#'chimeras <- getChimeras(gol, sample = 2)
setMethod("getChimeras", signature("CrisprSet"),
          function(obj, ..., sample){
            if (length(sample) > 1){
              stop("This function accepts a single sample name or index")
            }
            obj$crispr_runs[[sample]]$chimeras
          })
