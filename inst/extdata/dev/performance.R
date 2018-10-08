#performance.R
#

#Scenario
#function to split a string into codon
#input e.g. "GCATTCTTAGCATTT"
#OUTPUT:"GCA","TTC","TTA","GCA",TTT"
s="GCATTCTTAGCATTT"

substring(s,c(1,4),c(3,6))

codons_1<-function(x){
  first=nchar(x)
  first=seq(1,1,by =3)
  last=first+2
  return(substring(x,first,last))
}


unlist(strsplit(s,"(?<=...)",perl = TRUE))
#much slower
codons_2<-function(x){
  return(unlist(strsplit(x,"(?<=...)",perl = TRUE)))
}
library("stringi")
stri_extract_all_regex()

#very fast
codons_3<-function(x){
  return(unlist(stri_extract_all_regex(x,"(...)")))
}


#time the runtime of the function get to know which function if the best.

x= paste(sample(c("A","C","G","T"),99999,replace = TRUE),collapse="")


#==========
Sys.time()
as.integer(Sys.time())


start_time <- Sys.time()
output<-codons_3(x)
end_time <- Sys.time()
end_time - start_time

system.time({output<- codons_3(x)})
system.time({output<- codons_3(x)})


install.packages("microbenchmark")
install.packages("ggplot2")
library("microbenchmark")
microbenchmark(output=codons_1(x))
mb<- microbenchmark(output<- codons_3(x),times=1000)
library(ggplot2)
autoplot(mb)
