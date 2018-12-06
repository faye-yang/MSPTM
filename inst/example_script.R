#example scrip

if (! require(protViz, quietly=TRUE)) {
  install.packages("protViz")
  library(protViz)
}
if (! require(ggplot2, quietly=TRUE)) {
  install.packages("ggplot2")
  library(ggplot2)
}
if (! require(rTANDEM, quietly=TRUE)) {
  install.packages("rTANDEM")
  library(rTANDEM)
}
if (! require(colourpicker, quietly=TRUE)) {
  install.packages("colourpicker")
  library(colourpicker)
}


# example 1
#makrer ion
HexNAc_MarkerIons <- c(126.05495, 138.05495, 144.06552, 168.06552, 186.07608, 204.08665)
Glykan_MarkerIons <- c(109.02841, 127.03897, 145.04954, 163.06010, 325.11292)
ADP_Ribose <- c(136.0618, 250.0935, 348.0704, 428.0367)
#data
data(HexNAc)
#post translatonal modification
ptm0 <- cbind(AA="-", mono=0.0, avg=0.0, desc="unmodified", unimodAccID=NA)
ptm1 <- cbind(AA='N', mono=317.122300, avg=NA, desc="HexNAc", unimodAccID=2)
ptm2 <- cbind(AA='M', mono=147.035400, avg=NA, desc="Oxidation", unimodAccID=1)
m <- as.data.frame(rbind(ptm0, ptm1, ptm2))

#ptm<-hepler(data=HexNAc,modification=m,mZmarker_ions=HexNAc_MarkerIons)
intensity_plot(data=HexNAc,modification = m,mZmarker_ions=HexNAc_MarkerIons,search_engine="Mascot")



#example 2: rTandem search for mouse: there are 150810 proteins so it will at most have 150810 plots 
#To make the program faster you can choose to see only a subset of the plots

modification<-data.frame("type"=c("Carbamidomethyl","Oxidation"),
                         "monomass"=c(57.022, 16.0), "AA"=c("C","M"))

unzip("./inst/extdata/output_mouse.2018_12_04_19_57_17.t.xml.zip",exdir="./inst/extdata")
result.file <- "./inst/extdata/output_mouse.2018_12_04_19_57_17.t.xml"
uids<-c(12,2,731)
result <- GetResultsFromXML(result.file)
data<-tandem_get_data(result,modification,uids)
write.csv(data,paste('./inst/extdata/data_tandem.csv', sep=""),row.names=FALSE)
intensity_plot(data,modification,mZmarker_ions, search_engine="Tandem")


# [END]