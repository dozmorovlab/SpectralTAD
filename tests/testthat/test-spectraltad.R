context("test-spectraltad")

test_that("TAD calling works", {
  expect_equal(SpectralTAD(rao_chr2), 4)
})

