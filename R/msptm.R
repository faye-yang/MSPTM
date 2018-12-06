
#' \code{tandem_get_data} extract required data(modification, peptide sequence, mass 
#' to charge) from rtandem output file.
#'
#' transform the rtandem output file to a data variables that can be accepted by 
#' PTM_finderfunction.
#' @param result a output file from rtandem package that contain mass to charge,
#'  peptide sequence, protein information, modification
#' @param modification modification information that are investigated
#' @param uids proteins id that you want to investigate
#' @return A numeric vector.
#'
#' modification<-data.frame("type"=c("Carbamidomethyl","Oxidation")
#' ,"monomass"=c(57.022, 16.0), "AA"=c("C","M"))
#' tandem_get_data(result.file,modification)
#' @export
tandem_get_data<-function(result,modification,uids){
  
  #get PTM of interest
  #unique uid
  ptm_info <- data.frame(Uid=c(), Type = c(),At=c(),Sequence=c(), Pep_ids=c())
  unique_uids<-unique(uids)
  for(uid in unique_uids){
    #update dataframe
    ptm_info<-rbind(ptm_info,get_modified(uid,modification,result))
    
  }
  return(ptm_info)
  
}



#' \code{get_modified} extract required modification information (modification location,
#'  modified sequence, type of modification) from protein information
#'
#' @param uid prot.uid attribute from result.file's peptide slot
#' @param modification interested modifications
#' @param result.file full information of x!tandem search engine return
#' @return A dataframe that contain protein's uid, sequence, location of modification, 
#' type of modification.
#'  
#'  modification<-data.frame("type"=c("Carbamidomethyl","Oxidation"),
#'  "monomass"=c(57.022, 16.0), "AA"=c("C","M"))
#' get_modified(uid=4,modification)

get_modified<-function(uid,modification,result.file){
  #number of modification
  peptides <- GetPeptides(protein.uid=uid, results =result.file, expect =0.05 )
  check<-is.na(peptides$ptm.type)
  peptides<-peptides[!check]
  #all NA so do not analyze it
  if(all(check)){
   
    ptm <- data.frame(Uid=c(), Type = c(),At=c(),Sequence=c(), Pep_ids=c())
    return(ptm)
  }
  peptides <- peptides[peptides$ptm.type == modification$AA]
  pep_ids<-peptides$pep.id
  #initiate all attributes
  protein_id<-c()
  at<-c()
  type<-c() 
  sequence<-c()
  pep_ids<-c()
  num_modification <-length(peptides$ptm.at)
  for(i in 1:num_modification){
    #update empty lists
    if(is.null(at)){
      #update attributes
      protein_id<-c(protein_id,uid)
      index<- peptides$ptm.at[i]-peptides$start.position[i] +1
      at<-c(at,index)
      type<-c(type,peptides$ptm.type[i])
      sequence<-c(sequence,peptides$sequence[i])
      pep_ids<-c(pep_ids,peptides$pep.id[i])
    }
    else{
      index<- peptides$ptm.at[i]-peptides$start.position[i] +1 
      check_dup_loca<-(at==index)
      #one location can only have one modification.
      if(!any(check_dup_loca)){
        #update attributes
        protein_id<-c(protein_id,uid)
        #the modification index of sequence: modification index of 
        #whole sequence- its start index
        at<-c(at,index)
        type<-c(type,peptides$ptm.type[i])
        sequence<-c(sequence,peptides$sequence[i])
        pep_ids<-c(pep_ids,peptides$pep.id[i])
        
      }
    }
  }
  
  ptm<-data.frame(Uid=protein_id, Type = type,At=at,Sequence=sequence, Pep_ids=pep_ids)
  return(ptm)
  
}


# [END]