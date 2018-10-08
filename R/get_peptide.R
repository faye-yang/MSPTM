
#' \code{get_peptide} output the peptides information of the root of xml file.
#'A helper function of \code{tandem_get_data}
#'
#' @param root a root of xml file
#' @return A character vector.
#' @export
get_peptide <- function(root, index) {
    peptide = c()
    if (!is.null(root[[index]][[1]][[3]][[1]])) {
        peptide_raw = xmlToList(root[[index]][[1]][[3]][[1]])
        peptide_raw1 = as.list(strsplit(peptide_raw, "\n\t")[[1]])
        
        peptide_raw2 = c()
        for (pep in peptide_raw1) {
            temp = as.list(strsplit(pep, " ")[[1]])
            peptide_raw2 = c(peptide_raw2, temp)
            
        }
        size = length(peptide_raw2)
        peptide_raw2[size] = as.list(strsplit(peptide_raw2[[size]], "\n")[[1]])[1]
        peptide = c(peptide, peptide_raw2)
    }
    
    return(peptide)
}


