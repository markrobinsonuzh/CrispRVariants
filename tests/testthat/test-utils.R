context("Counting mutations")

local({
    alns <- GenomicAlignments::GAlignments(seqnames=rep("1", 5), 
                                       pos=as.integer(rep(1,5)), 
                    cigar=c("10M", "5M1D5M", "4M1D6M4D2M",
                            "1M2D3M4I", "2M3I2M3I2M"), 
                    strand = S4Vectors::Rle(factor(rep("+", 5),
                                         levels = c("+","-","*"))))

    test_that("countDeletions returns the expected counts", {
      expect_equal(countDeletions(alns), 1)
      expect_equal(countDeletions(alns, multi.del = TRUE), 2)
      expect_equal(countDeletions(alns, del.and.ins = TRUE), 2)
      expect_equal(countDeletions(alns, multi.del = TRUE, del.and.ins = TRUE), 3) 
      expect_equal(countDeletions(alns, del.ops=c("N")), 0) 
    })

    test_that("countInsertions returns the expected counts", {
      expect_equal(countInsertions(alns), 0)
      expect_equal(countInsertions(alns, multi.ins = TRUE), 1)
      expect_equal(countInsertions(alns, ins.and.del = TRUE), 1)
      expect_equal(countInsertions(alns, multi.ins = TRUE, ins.and.del = TRUE), 2) 
    })

    test_that("countIndels returns the expected counts", {
      expect_equal(countIndels(alns), 4)
    })

    test_that("indelPercent returns the expected value", {
      expect_equal(indelPercent(alns), (4/5)*100)
    })
})