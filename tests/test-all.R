library(testthat)
library("crispRvariants")
test_package("crispRvariants")

# Test CrisprSet: filterUniqueLowQual - not sure if it's removing 
# all sequences in a row when only one matches delete criteria

# Test CrisprSet: countVariantAlleles - not sure that it actually is 
# excluding SNVs

#target <- GRanges("chr1", IRanges(20, 30))
#cigars <- c("5M6D10M","2M10I3M6D10M","5M6D5M","2M5I3M6D5M")
#starts <- c(13,13,15,15)
