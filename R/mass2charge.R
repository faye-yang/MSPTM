
#' get mass to charge value from thexml root.
#' @param rtandem_oput a output file from rtandem package
#' @return A numeric vector.
#' @export
get_mass2charge <- function(root, index) {
  mz_list = c()
  if (!is.null(root[[index]][[1]][[3]][[1]])) {
    mz_raw = xmlToList(root[[index]][[3]][[2]])
    if (typeof(mz_raw) == "list") {
      attri = mz_raw$Xdata$.attrs
      # only for the section that is to find mass to charge
      grep_list = grepl("MASSTOCHARGERATI", attri)
      counter = 0
      for (boo in grep_list) {
        
        if (boo == TRUE) {
          counter = 1
        }
        if (counter == 1) {
          
          mz_split = mz_raw$Xdata$values$text
          mz = as.list(strsplit(mz_split, "\n")[[1]])
          for (i in mz) {
            # print(i)
            if (i != "") {
              temp = as.list(strsplit(i, " ")[[1]])
              mz_list = c(mz_list, temp)
            }
          }
        }
      }
      index = index + 1
    }
  }
  return(mz_list)
  
}
