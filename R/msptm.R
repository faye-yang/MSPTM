#install.packages('devtools')
library(roxygen2)

# input: a list of mz, modification marker


#' transform the rtandem output file to a data variables that can be accepted by PTM_finderfunction.
#' @param rtandem_oput a output file from rtandem package
#' @return A numeric vector.
#' @export
tandem_get_data<-function(rtandem_oput){
  if( file_ext(rtandem_oput) !="xml" ) stop('input file is not a xml file')
  rtandem_data<- xmlParse(file = rtandem_oput)
  root=xmlRoot(rtandem_data)
  size=xmlSize(root)
  if( size ==0 ) stop('input file has zero information to analyze')
  datalist=list()
  for(i in 1:size){
    
    #get mass to charge
    mz_list=get_mass2charge(root,i)
    
    #get peptide
    peptide_list=get_peptide(root,i)
    
    if(!(is.null(peptide_list) & is.null(mz_list))){
      data_tandem=list(peptide=peptide_list,mz=mz_list)
      datalist=c(datalist, data_tandem)
    
    }
  }
  return(datalist)
  
}



