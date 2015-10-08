context("Narrowing alignments to a target region")

# Setup data
cigars <- c("5M5D20M","10M10D20M","5M10D20M","5M20D10M","1M3I2M3I3M10D10M",
            "10M5I30M", "15M2D20M", "5M", "12M5I15M", "12M9I10M", "3M6D15M",
            "4S3M1I6M4D8M")
wdths <- GenomicAlignments::cigarWidthAlongPairwiseSpace(cigars)
gr <- GenomicRanges::GRanges("chr1", IRanges(start = 1, width = wdths), 
                             cigar = cigars, strand = "+")
galns <- as(gr, "GAlignments")
target <- GenomicRanges::GRanges("chr1", IRanges(11, 20))
expected <- c("10M", "10D", "10D5M", "20D", "10D4M", "10M", 
              "5M2D3M", "2M5I8M","2M9I8M","10M", "4D7M")
narrowed <- narrowAlignments(galns, target)


# Test 1
test_that("narrowAlignments correctly narrows to the nearest operation",{
  expect_equal(cigar(narrowed), expected)
})


# Test 2
quals <- Biostrings::DNAStringSet(sapply(wdths, function(x){
  paste0(rep("H",x), collapse = "")})) 
mcols(galns)$qual <- quals

test_that("narrowAlignments with qual but no seq produces an error",{
  expect_error(narrowAlignments(galns, target))
})


# Test 3
seqs <- Biostrings::DNAStringSet(sapply(wdths, function(x) paste0(rep("A",x), collapse = ""))) 
mcols(galns)$seq <- seqs
narrowed <- narrowAlignments(galns, target)

test_that("narrowAlignments narrows quality scores",{
  expect_equal(width(mcols(narrowed)$qual), width(mcols(narrowed)$seq))
})


# Test 4
mcols(galns)$dat <- seq_along(galns)
test_that("narrowAlignments preserves metadata",{
  expect_equal(length(mcols(narrowAlignments(galns, target))), 3)
})
