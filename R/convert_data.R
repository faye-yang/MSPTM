
#' refereence: perl script comes from package protViz get rdatafile from a mascot output file
#'
#' @param iput_file input file that is received from mascot
#' @param oput_file a file that is readable by the protviz package
#' @export
mascot_data<- function(iput_file,oput_file){
  #mascot data to rdata file 
  get_data <- paste("perl", "mascotDat2RData.pl", "-d="+iput_file, "-m="+oput_file)
  system(get_data)
  return(oput_file)
}

# [END]