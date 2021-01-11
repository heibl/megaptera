## Last update 2020-0223
library(testthat)
library(megaptera)

## Test strip.infraspec()
test_that("return genus and epitheton", {
  expect_identical(strip.infraspec("Pinus mugo subsp. x rotundata voucher 06_0526"), "Pinus mugo")
          })
## Test indet.strings()
# indet <- indet.strings(hybrids = FALSE, collapse = TRUE)
# test_that("undetermined names are identified correctly", {
#   expect_length(grep(indet, "Pinus mugo"), 0)
#   expect_length(grep(indet, "Pinus_mugo"), 0)
#   expect_length(grep(indet, "Pinus"), 1)
#   expect_length(grep(indet, "Onthophagus nr. babirussa"), 1) # nr = "near"
#   expect_length(grep(indet, "Amanita_sp"), 1)
#   expect_length(grep(indet, "Amanita_sp."), 1)
#   expect_length(grep(indet, "Amanita_sp_xxx"), 1)
#   expect_length(grep(indet, "Amanita_sp-53"), 1)
#   expect_length(grep(indet, "Hydropsyche_n.sp._2006031401"), 1)
#   expect_length(grep(indet, "Amylosporus_sp._'succulentus'"), 1)
#   expect_length(grep(indet, "Limenitis_hybrid_form_'rubidus'"), 1)
#   expect_length(grep(indet, "Juniperus_hybrid"), 1)
#   expect_length(grep(indet, "Juniperus_hybrid_sp._LO-2009"), 1)
#   expect_length(grep(indet, "Juniperus_hybridus"), 0)
# })


## Test is.Linnean()
test_that("identify valid scientific binomials", {
  expect_true(is.Linnean("Pinus mugo"))
  expect_true(is.Linnean("Pinus mugo subsp. x rotundata"))
  # expect_false(is.Linnean("Pinus mugo subsp. x rotundata voucher 06_0526")) ## IMPLEMENT ME!
  expect_true(is.Linnean("Myotis riparius PS3"))
  expect_false(is.Linnean("Vespertilionidae"))
})

## Test splitGiTaxon()
test_that("correctly split taxon name and identifier", {
  ## BOLD style identifier
  expect_identical(splitGiTaxon("Strix nigrolineata USNMA071_10"), 
                   list(taxon = "Strix nigrolineata", gi = "USNMA071_10"))
  expect_identical(splitGiTaxon("Strix_nigrolineata_USNMA071_10"), 
                   list(taxon = "Strix_nigrolineata", gi = "USNMA071_10"))
  expect_identical(splitGiTaxon("Strix nigrolineata_USNMA071_10"), 
                   list(taxon = "Strix nigrolineata", gi = "USNMA071_10"))
  expect_identical(splitGiTaxon("Strix_nigrolineata USNMA071_10"), 
                   list(taxon = "Strix_nigrolineata", gi = "USNMA071_10"))
  ## NCBI style identifier
  expect_identical(splitGiTaxon("Bufo bufo AY840247"), 
                   list(taxon = "Bufo bufo", gi = "AY840247"))
  ## epitheton with hyphan
  expect_identical(splitGiTaxon("Bufo bufo-bufo AY840247"), 
                   list(taxon = "Bufo bufo-bufo", gi = "AY840247"))
  ## subspecies
  expect_identical(splitGiTaxon("Strix uralensis macrourus AY840247"), 
                   list(taxon = "Strix uralensis macrourus", gi = "AY840247"))
})

