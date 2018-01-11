library(testthat)
library(megaptera)

test_that("return genus and epitheton", {
  expect_identical(strip.infraspec("Pinus mugo subsp. x rotundata voucher 06_0526"), "Pinus mugo")
          })


test_that("identify valid scientific binomials", {
  expect_true(is.Linnean("Pinus mugo"))
  expect_true(is.Linnean("Pinus mugo subsp. x rotundata"))
  # expect_false(is.Linnean("Pinus mugo subsp. x rotundata voucher 06_0526")) ## IMPLEMENT ME!
  expect_true(is.Linnean("Myotis riparius PS3"))
  expect_false(is.Linnean("Vespertilionidae"))
})

