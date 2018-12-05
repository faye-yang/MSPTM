#msptm.R
devtools::install_github("MangoTheCat/visualTest")
library(visualTest)
context("msptm")

# ==== BEGIN SETUP AND PREPARE =================================================

#get the package to the right address
#alt(option) + cursor(shubiao) can select multiple line or eddit multiple line at the same time



#
# ==== END SETUP AND PREPARE ===================================================

test_that("corrupt input generates errors from helper function get_peptide()",  {
  expect_error(intensity_plot(), "argument \"search_engine\" is missing, with no default")
  expect_error(tandem_get_data(),"argument \"uids\" is missing, with no default")
})



library(rTANDEM)



# [END]
