library(ggplot2)

#' \code{dist_colour_plot} output graph of .
#'A helper function of \code{tandem_get_data}
#'
#' @param rdata mass spetrometry information for the peptide
#' @param modification contain modification information , intensity of ion, amino acide that is modified
#' @param mZmarkerIons maker ion 
#' @export
intensity_plot<- function(data,modification,mZmarker_ions){
  ptm<-PTM_MarkerFinder(
    data,
    modification$mono,
    modification$desc,
    mZmarker_ions,
    minNumberIons=2,
    itol_ppm=10,
    minMarkerIntensityRatio=5,
    PEAKPLOT =TRUE
  )
  dev.off()
  #print(ptm)
  
  #find the peptides that are modified
  
  for (peptide_inf in data){
    peptide <- peptide_inf$peptideSequence
    num_aa<-length(peptide)
    string2char<-strsplit(peptide, "")
    #initialize the intensity for each aa
    num_aa<-length(string2char[[1]])
    intensity=rep(0, num_aa)
    for(aa_index in c(1:num_aa) ){
      num_modi<-length(modification$desc)
      aa<-string2char[[1]][aa_index]
      for(index in c(1:num_modi)){
        if(aa==as.character(modification$AA[index]) & (aa!="-")){
            #print("peptide:",peptide,"ptm object peptide:",ptm$peptideSequence)
            intensity[aa_index]<-ptm$markerIonIntensity[ptm$peptideSequence==peptide][num_modi]
            print(peptide)
            print(ptm$markerIonIntensity)
            
          }
        }

    }
    #plot:only the position has modification will have marker intensity
    if( peptide %in% ptm$peptideSequence ){
    peptide_frame<-data.frame(pep_sequence=string2char[[1]],markerIonIntensity=intensity)
    print(peptide_frame)
    p<- ggplot(ptm, aes(x=ptm$pep_sequence,y=ptm$markerIonIntensity))+ 
        #geom_bar(stat="identity", width = 0.5, fill="tomato2") + 
        labs(title=paste("Post translational modification intensity 
                         for peptide", peptide, sep=" "))
    print(p)
    }
    
   
  }
}



# [END]