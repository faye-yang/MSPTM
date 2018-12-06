#' \code{intensity_plot} output graph of intensity of amino acid.
#'
#' @param data mass spetrometry information for the peptide
#' @param modification contain modification information , intensity of ion, amino acide that is modified
#' @param mZmarker_ions maker ion
#' @param search_engine can be Mascot or Tandem
#' @examples
#' modification<-data.frame("type"=c("Carbamidomethyl","Oxidation"),
#' "monomass"=c(57.022, 16.0), "AA"=c("C","M"))
#' result.file <- "/Users/yufei/Desktop/2018fall/BCB410/MSPTM/inst/extdata/output_mouse.2018_12_04_19_57_17.t.xml"
#' uids<-c(12,2,731)
#' library(rTANDEM)
#' result <- GetResultsFromXML(result.file)
#' data<-tandem_get_data(result,modification,uids)
#' intensity_plot(data,modification,mZmarker_ions, search_engine="Tandem")

#' @export
intensity_plot<- function(data,modification,mZmarker_ions, search_engine){
  if(search_engine == "Mascot"){
    print("Mascot")
    ptm<-PTM_MF(
      data,
      modification$mono,
      modification$desc,
      mZmarker_ions,
      minNumberIons=2,
      itol_ppm=10,
      minMarkerIntensityRatio=5,
      PEAKPLOT =TRUE
    )
    
    #find the peptides that are modified
    num_modified<-length(ptm$modification)
    for(i in 1:num_modified){
      #get modification indicators
      modified<- as.character(ptm$modification[i])
      sequence<-as.character(ptm$peptideSequence[[i]])
      print(sequence)
      if(!sequence=="NA"){
        # if indicator is 1. index of indicator is the index of modified aa
        modified_index<-as.numeric(unlist(strsplit(modified,"")))
        end<-length(modified_index)
        modified_index<-modified_index[c(-1,-end)]
        seq<-as.character(unlist(strsplit(sequence,"")))
        #create data frame for sequence and intensity
        mascot_data<-data.frame(modified_index,seq)
      
        #get intensity
        #update intensity with correct index
        avg_intensity<-mean(ptm$markerIonIntensity)
        p<-ggplot(data=mascot_data, aes(x=seq, y=modified_index*avg_intensity)) +
          geom_bar(stat="identity", fill="steelblue")+
          theme_minimal()
        #display plot
        print(p)
      }
    }
  }
  
  if(search_engine =="Tandem"){
    print("plot Tandem ptm graph")
    rows<-nrow(data)
    for(i in 1:rows){
      sequence<-as.character(data$Sequence[[i]])
      seq<-as.character(unlist(strsplit(sequence,"")))
      size<-length(seq)
      indicators<-c(rep(0, size))
      indicators[data$At[i]]<-1
      #create dataframe for PTM seq and indicators
      data_PTM<-data.frame(indicators,seq)
      #get protein and peptide id
      pep_id<-data$Pep_ids[[i]]
      uid<-data$Uid[[i]]
      #save data as .csv file for rshiny to plot 
      write.csv(data_PTM,paste('./inst/extdata/data_PTM',as.character(i),
                              '.csv', sep=""),row.names=FALSE)
      #plot the graph
      title<-paste("Protein sequence:",uid," with peptide id: ",pep_id,"vs PTM")
      p<- ggplot(data=data_PTM, aes(x=seq, y=indicators)) + 
        geom_bar(stat="identity", fill="steelblue")+
        theme_minimal() +
        labs(title=title, y="Post-translational Modification", x="Protein sequence")
      #print out the plot
      print(p)
 
    }
  }


   
}

# [END]