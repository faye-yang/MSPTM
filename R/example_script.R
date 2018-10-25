#example scrip

install.packages("protViz")
library(protViz)

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
intensity_plot(data=HexNAc,modification = m,mZmarker_ions=HexNAc_MarkerIons)
