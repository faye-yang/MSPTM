
#' \code{dist_colour_plot} output graph of .
#'A helper function of \code{tandem_get_data}
#'
#' @param rdata mass spetrometry information for the peptide
#' @param modificationName name of intensity
#' @param modification intensity of ion
#' @param mZmarkerIons maker ion 
#' @export
dist_colour_plot<- function(rdata,modificationName,modification,mZmarkerIons){
  ptm<-PTM_MarkerFinder(
    rdata,
    modification,
    modificationName,
    mZmarkerIons,
    minNumberIons=2,
    itol_ppm=10,
    minMarkerIntensityRatio=5,
    mgfFilename=-1,
    PEAKPLOT =FALSE
  )
  num_color=1
  for (i in modificationName) {
    mod_list=list()
    if (i !=":unmodified"){
      if(!is.element(i, mod_list)){
        mod_list[[num_color]]=i
        num_color=num_color+1
      }
    }
  }
  colour<-distinctColorPalette(num_color)
  x <- modification
  y <- mZmarkerIons
  plot(x,y, main="scatter plot of post translational modification and its intensity", 
       col=colour, pch=16)
  dev.off()
  
}

