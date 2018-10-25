# MSPTM


## R Package MSPTM

This is a RStudio project for R packages,that is use to visualize the quantification of
post translational modifcation of protein in mass spectrometry.
-----------------------------------------------

Note: you can't push empty directories to your repository. Make sure youu keep
at least one file in every directory that you want to keep during development.

-----------------------------------------------

Some useful keyboard shortcuts for package authoring:

* Build and Reload Package:  `Cmd + Shift + B`
* Update Documentation:      `Cmd + Shift + D` or `devtools::document()`
* Test Package:              `Cmd + Shift + T`
* Check Package:             `Cmd + Shift + E` or `devtools::check()`

-----------------------------------------------


Load the package (outside of this project) with:
`devtools::install_github("faye-yang/MSPTM")`

background information:
X!Tandem and mascot are the two most popular database enginge to search PTM
rTandem interfaces the X!Tandem protein identification algorithm. mascot is a database search engine to search PTM and provide many other mass spectrometry information to perform PTM identification.


convert_data.R: reference the perl file from Protviz package and convert mascot data to rdata
get_peptide.R: It contains an R function that extract peptide sequence from the rTanem engine output
mass2charge.R:It contains an R function that extract mz from rTandem engine output
msptm.R : It contains an R function that create a data frame for the rTandem engine input for plotting input( incomplete)
plot.R: It contains an R function that can plot the peptide modification plot.

