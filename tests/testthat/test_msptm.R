#msptm.R
devtools::install_github("MangoTheCat/visualTest")
library(visualTest)
context("msptm")

# ==== BEGIN SETUP AND PREPARE =================================================
load(system.file("extdata/testdata",".Rdata",package="msptm")) 
#get the package to the right address
#alt(option) + cursor(shubiao) can select multiple line or eddit multiple line at the same time

#tmp <- c(  1.0
#         , 1.778279410038922980775
#         , 3.162277660168379522787
#         , 5.623413251903492060535
#         , 10.0)

#
# ==== END SETUP AND PREPARE ===================================================

test_that("corrupt input generates errors from helper function get_peptide()",  {
  expect_error(get_peptide(), "argument \"x\" is missing, with no default")
  expect_error(get_peptide(root,"i"),"input type is not correct")
})
test_that("corrupt input generates errors from function tandem_get_data()",  {
  expect_error(tandem_get_data(), "argument \"x\" is missing, with no default")
  expect_error(tandem_get_data("output1"),"input format incorrect")
})


rtandem_data<- xmlParse(file = "/Users/yufei/Desktop/2018fall/BCB410/MSPTM/output1.xml")
root=xmlRoot(rtandem_data)
tmp=NULL
rdata=tandem_get_data("/Users/yufei/Desktop/2018fall/BCB410/MSPTM/output1.xml")
test_that("a sample input prouces the expected output from helper function get_peptide()",  {
  expect_equal(get_peptide(root,12), tmp)
  
})


test_that("a sample input prouces the expected output",  {
  expect_equal(lseq(1,10, length.out = 5), tmp)
  
})



# ==== BEGIN TEARDOWN AND RESTORE ==============================================
# Remove every persitent construct that the test has created, except for
# stuff in tempdir().
#
rm()

# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
