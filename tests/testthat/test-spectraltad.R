context("test-spectraltad")

test_that("TAD calling works", {
  data("rao_chr20_25_rep")
  expect_equal(nrow(SpectralTAD(rao_chr20_25_rep, "chr20")$Level_1), 165)
})

