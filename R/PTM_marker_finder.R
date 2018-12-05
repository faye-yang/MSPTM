

# reference: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/PTM_MarkerFinder.R $
# $Id: PTM_MarkerFinder.R 6318 2014-04-02 13:02:51Z cpanse $
# modified from protviz without the error and adjust the margin


#' \code{PTM_MF} output graph of .
#'A helper function of intensity_plot
#'
#' @param data mass spetrometry information for the peptide
#' @param modification contain modification information , intensity of ion, amino acide that is modified
#' @param modi_name name of modification
#' @param mZmarkerIons maker ion 
#' @param itol_ppm ppm
#' @param minMarkerIntensityRatio minimum ratio for marker ion intensity 
#' @param mgfFilename mgf file name indicator
#' @param PEAKPLOT peak plot for MTP 
#' @param minNumberIons minimum number or marker ion


PTM_MF<-function(data, 
                   modification, 
                   modi_name, 
                   mZmarkerIons, 
                   minNumberIons=2, 
                   itol_ppm=10, 
                   minMarkerIntensityRatio=5,
                   mgfFilename=-1, 
                   PEAKPLOT=TRUE){
  
  if(!PEAKPLOT){
    .PTM_MarkerFinder_no_peakplot(data, 
                                  modification,
                                  modi_name,
                                  mZmarkerIons, 
                                  minNumberIons, 
                                  itol_ppm,
                                  minMarkerIntensityRatio)
    return (TRUE)
  }
  
  
  query.idx<-1:length(data)
  query.to.scan<-as.integer(as.character(lapply(data, function(x){
    if (length(x$scans) == 1){
      return(x$scans)
    }else{return (x$scans[1])}
  })))
  
  scan.to.query<-rep(-1, max(query.to.scan, na.rm=TRUE))
  scan.to.query[query.to.scan[query.idx]]<-query.idx
  
  
  if (mgfFilename != -1){
    unlink(mgfFilename)
    .PTM_MarkerFinder_writeMGF_Header(mgfFilename)
  }
  
  rr<-numeric()
  
  for (i in 2:length(data)){
    
    idx <- findNN_(mZmarkerIons, data[[i]]$mZ)
    
    # ppm.itol.cutoff
    ppm.error <- 1e-06 * itol_ppm * data[[i]]$mZ[idx]
    b<-(abs(mZmarkerIons-data[[i]]$mZ[idx]) < ppm.error)
    
    sum.mZmarkerIons.intensity <- sum(data[[i]]$intensity[idx[b]])
    sum.intensity <- sum(data[[i]]$intensity) - sum.mZmarkerIons.intensity
    percent.mZmarkerIons <- round(100 * (sum.mZmarkerIons.intensity / (sum.mZmarkerIons.intensity + sum.intensity)),1)
    
    if ((length(data[[i]]$mZ[idx[b]]) >= minNumberIons) 
        & percent.mZmarkerIons > minMarkerIntensityRatio ){
      
      if (mgfFilename != -1){
        
        .PTM_MarkerFinder_writeMGF(data[[i]], mgfFilename, 
                                   pattern=paste(data[[i]]$mZ[idx[b]], data[[i]]$intensity[idx[b]]))
      }
      
      r <- cbind(scans=data[[i]]$scans, 
                 mZ=data[[i]]$mZ[idx[b]],
                 markerIonMZ=mZmarkerIons[b],
                 markerIonIntensity=data[[i]]$intensity[idx[b]], 
                 markerIonMzError=mZmarkerIons[b]-data[[i]]$mZ[idx[b]],
                 markerIonPpmError=1e+06 * (mZmarkerIons[b]-data[[i]]$mZ[idx[b]])/data[[i]]$mZ[idx[b]],
                 query=i,
                 pepmass=data[[i]]$pepmass,
                 peptideSequence=data[[i]]$peptideSequence,
                 modification=data[[i]]$modification
      )
      
      
      rr<-rbind(rr,r)
      
      
      
      ####### P E A K P L O T ########################################################
      if (!is.na(data[[i]]$mascotScore) ){
        fi<-fragmentIon(sequence=data[[i]]$peptideSequence, 
                        FUN=defaultIon,
                        modified=substr(data[[i]]$modification, 2, nchar(data[[i]]$modification)-1),
                        modification=modification)
        
        fi.by<-as.data.frame(cbind(b=fi[[1]]$b, y=fi[[1]]$y))
        par(mar=c(1,1,1,1))
        peakplot(data[[i]]$peptideSequence, 
                 spec=data[[i]], fi=fi.by, ion.axes=FALSE,  
                 main=paste("scantype: HCD / peptide sequence: ", 
                            .PTM_MarkerFinder_(data[[i]]$peptideSequence, data[[i]]$modification, modi_name), sep=''),
                 xlim=c(0,max(data[[i]]$mZ)))
       
      }else{
        par(cex=0.75,mar=c(1,1,1,1))
        plot(data[[i]]$mZ, data[[i]]$intensity, type='h',
             xlab='m/z',ylab='Intensity',
             main=paste("scantype: HCD"),
             xlim=c(0, max(data[[i]]$mZ)),
             sub=paste(i, data[[i]]$title))
      }
      
      points(data[[i]]$mZ[idx[b]],
             data[[i]]$intensity[idx[b]], 
             pch=22,
             col='#6CC417AA', 
             bg='#6CC417AA', 
             cex=1)
      
      legend("topright", paste(c('protein', 'm/z', 'charge', 'ionScore', 'query#','xlab','ylab'), 
                               c(data[[i]]$proteinInformation, 
                                 data[[i]]$pepmass, data[[i]]$charge, data[[i]]$mascotScore,i,'m/z','intensity')), cex=0.75)
      
      
      ################################################################################
      ####### P E A K P L O T ########################################################
      j<-1;
      PLOTFLAG<-FALSE
      while (j<5){ 
        #k<-i+j
        k <- scan.to.query[query.to.scan[i] + j]
        
        if (k < 0){break;}
        
        if (!is.na(data[[k]]$mascotScore) & abs(data[[k]]$pepmass - data[[i]]$pepmass) < 1){
          
          
          if (mgfFilename != -1){
            .PTM_MarkerFinder_writeMGF(data[[k]], mgfFilename)
          }
          
          fi<-fragmentIon(sequence=data[[k]]$peptideSequence, 
                          FUN=defaultIon,
                          modified=substr(data[[k]]$modification, 2, nchar(data[[k]]$modification)-1),
                          modification=modification)
          
          fi.cyz<-as.data.frame(cbind(y=fi[[1]]$y, c=fi[[1]]$c, z=fi[[1]]$z))
          par(mar=c(1,1,1,1))
          p<-peakplot(data[[k]]$peptideSequence, spec=data[[k]], 
                      main=paste("scantype: ETD / peptide sequence: ",
                                 .PTM_MarkerFinder_(data[[k]]$peptideSequence, data[[k]]$modification, modi_name), sep=''),
                      fi=fi.cyz,
                      itol=0.6,
                      xlim=c(0, max(data[[i]]$mZ)),
                      ion.axes=FALSE)
          PLOTFLAG<-TRUE
          
          
          legend("topleft", paste(c('protein', 'm/z', 'charge', 'ionScore', 'query#','xlab','ylab'), 
                                  c(data[[k]]$proteinInformation, 
                                    data[[k]]$pepmass, data[[k]]$charge, data[[k]]$mascotScore, k, 'm/z','intensity')), cex=0.75)
          break;
        }
        j <- j+1
      }
      
      if (PLOTFLAG == FALSE){
        par(mar=c(1,1,1,1))
        
        plot(0,0,type='n', xlab='',ylab='', axes=FALSE)
        text(0,0,"no peptide assignment", cex=3, col="lightgrey")
      }
      ################################################################################
      
      ################################################################################
      
      par(mar=c(1,1,1,1))
      plot(mZmarkerIons, (mZmarkerIons-data[[i]]$mZ[idx]),
      plot(mZmarkerIons[b], 1e+06 * (mZmarkerIons[b]-data[[i]]$mZ[idx[b]])/data[[i]]$mZ[idx[b]],
           axes=TRUE,type='p', xlab='m/z',ylab='ppm error',log=''), ylim=c(-max(ppm.error),max(ppm.error)))
      axis(3, mZmarkerIons[b], round(mZmarkerIons[b],1),las=2)
      abline(h=0.0, col='grey')
      legend("topright", paste(c('xlab','ylab'), 
                               c('m/z for maker Ion','ppm error')), cex=0.75)
      
      
      
    
    }
  }
  
  # TODO(cp): make a S3 class
  rr <- as.data.frame(rr)
  
  rr$markerIonIntensity <- as.numeric(levels(rr$markerIonIntensity))[rr$markerIonIntensity]
  rr$mZ <- as.numeric(levels(rr$mZ))[rr$mZ]
  rr$pepmass <- as.numeric(levels(rr$pepmass))[rr$pepmass]
  rr$markerIonMZ <- as.numeric(levels(rr$markerIonMZ))[rr$markerIonMZ]
  rr$markerIonMzError <- as.numeric(levels(rr$markerIonMzError))[rr$markerIonMzError]
  
  return(rr)
}

#' \code{.PTM_MarkerFinder_writeMGF_Header} a helper function to write the mgf header  .
#'
#' @param mgfFilename file name of the MGF 
#.PTM_MarkerFinder_writeMGF_Header(mgfFilename)
.PTM_MarkerFinder_writeMGF_Header <- function(mgfFilename){
  FILE <- file(mgfFilename, "a")
  writeLines(as.character(paste("# protViz packageVersion ", packageVersion('protViz'))), FILE)
  writeLines(paste("# ",date(),sep=''), FILE)
  close(FILE)
}
#' \code{.PTM_MarkerFinder_writeMGF} a helper function to find marker ion
#'
#' @param spec spectrum of peptides 
#' @param mgfFilename file name of the MGF 
#' @param pattern indicator for pattern plot
#' 

.PTM_MarkerFinder_writeMGF <- function(spec, mgfFilename, pattern=-1){
  FILE <- file(mgfFilename, "a")
  if (length(pattern)>0){
    writeLines(as.character(paste("# PATTERN ", pattern, sep=' ')), FILE)
  }
  writeLines("BEGIN IONS", FILE)
  writeLines(paste("TITLE=",as.character(spec$title),sep=''), FILE)
  writeLines(paste("PEPMASS=",as.character(spec$pepmass),sep=''), FILE)
  writeLines(paste("CHARGE=",as.character(spec$charge),'+',sep=''), FILE)
  writeLines(paste("SCANS=",as.character(spec$scans),sep=''), FILE)
  writeLines(paste("RTINSECONDS=",as.character(spec$rtinseconds),sep=''), FILE)
  writeLines(as.character(paste(spec$mZ, spec$intensity, sep=' ')), FILE)
  writeLines("END IONS", FILE)                              
  writeLines("", FILE)                              
  close(FILE)
}

#' \code{.PTM_MarkerFinder_} a helper function if peakplot option is FALSE .
#' @param modification contain modification information , intensity of ion, amino acide that is modified
#' @param modi_name name of modification
#' @param peptideSequence peptide sequence 

.PTM_MarkerFinder_ <- function(peptideSequence, modification, modi_name){
  
  result <- as.character()
  
  n <- nchar(peptideSequence)
  seq <- strsplit(peptideSequence,'')[[1]]
  mod <- as.numeric((strsplit(modification,'')[[1]]))
  
  for (i in 1:n){
    result<-paste(result, seq[i], sep='')
    if (mod[i+1] > 0){
      result<-paste(result,'[', modi_name[mod[i+1]+1], ']', sep='')
    }
  }
  return(result)
}
#' \code{.PTM_MarkerFinder_no_peakplot} a helper function 
#' for no peak plot option
#' @param data mass spetrometry information for the peptide
#' @param modification contain modification information , intensity of ion, amino acide that is modified
#' @param modi_name name of modification
#' @param mZmarkerIons maker ion 
#' @param minNumberIons minimum number of marker ion
#' @param itol_ppm ppm
#' @param minMarkerIntensityRatio minimum ratio for marker ion intensity 


.PTM_MarkerFinder_no_peakplot<-function(data, 
                                        modification,
                                        modi_name,
                                        mZmarkerIons, 
                                        minNumberIons=2, 
                                        itol_ppm=10, 
                                        minMarkerIntensityRatio=5){
  
  lapply(data, function(x){
    idx <- findNN_(mZmarkerIons, x$mZ)
    ppm.error <- 1e-06 * itol_ppm * x$mZ[idx]
    b<-(abs(mZmarkerIons-x$mZ[idx]) < ppm.error)
    
    sum.mZmarkerIons.intensity <- sum(x$intensity[idx[b]])
    sum.intensity <- sum(x$intensity) - sum.mZmarkerIons.intensity
    percent.mZmarkerIons <- round(100 * (sum.mZmarkerIons.intensity / (sum.mZmarkerIons.intensity + sum.intensity)),1)
    
    if ((length(x$mZ[idx[b]]) >= minNumberIons) & percent.mZmarkerIons > minMarkerIntensityRatio ){
      
      par(mar=c(1,1,1,1))
      plot(x$mZ, x$intensity, type='h', xlab='m/z',ylab='Intensity', xlim=c(0, max(x$mZ)), sub=paste(x$title))
      points(x$mZ[idx[b]], x$intensity[idx[b]], pch=22, col='#6CC417AA', bg='#6CC417AA', cex=1)
      legend("topright", paste(c('query', 'pepmass', 'charge','xlab','ylab'), c(x$id, x$pepmass, x$charge,'m/z','intensity')), cex=0.75)
      
      ################################################################################
      
      #plot(mZmarkerIons, (mZmarkerIons-data[[i]]$mZ[idx]),
      par(mar=c(1,1,1,1))
      plot(mZmarkerIons[b], 1e+06 * (mZmarkerIons[b] - x$mZ[idx[b]]) / x$mZ[idx[b]],
           axes=TRUE, type='p', xlab='m/z', ylab='ppm error',log='')# ylim=c(-max(ppm.error),max(ppm.error)))
      axis(3, mZmarkerIons[b], round(mZmarkerIons[b],1),las=2)
      abline(h=0.0, col='grey')
      
     
      #########  Pattern plot############
      par(mar=c(1,1,1,1))
      plot(data[[i]]$mZ[idx[b]],
           data[[i]]$intensity[idx[b]],
           col='#6CC417AA', 
           lwd=2, 
           type='h',
           xlab='m/z',
           ylab='Intensity',
           axes=TRUE)
      axis(3,data[[i]]$mZ[idx[b]],round(data[[i]]$mZ[idx[b]],1), las=2)
      legend("topright", paste(c('xlab','ylab'), c(,'m/z for maker ion','intensity')), cex=0.75)
    }
    
  }
  )
}


# [END]
