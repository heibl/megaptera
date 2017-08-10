library(testthat)
library(megaptera)

test_that("return genus and epitheton", {
  expect_identical(strip.infraspec("Pinus mugo subsp. x rotundata voucher 06_0526"), "Pinus mugo")
          })

