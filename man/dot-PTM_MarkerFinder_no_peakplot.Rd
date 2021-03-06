% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PTM_marker_finder.R
\name{.PTM_MarkerFinder_no_peakplot}
\alias{.PTM_MarkerFinder_no_peakplot}
\title{\code{.PTM_MarkerFinder_no_peakplot} a helper function 
for no peak plot option}
\usage{
.PTM_MarkerFinder_no_peakplot(data, modification, modi_name, mZmarkerIons,
  minNumberIons = 2, itol_ppm = 10, minMarkerIntensityRatio = 5)
}
\arguments{
\item{data}{mass spetrometry information for the peptide}

\item{modification}{contain modification information , intensity of ion, amino acide that is modified}

\item{modi_name}{name of modification}

\item{mZmarkerIons}{maker ion}

\item{minNumberIons}{minimum number of marker ion}

\item{itol_ppm}{ppm}

\item{minMarkerIntensityRatio}{minimum ratio for marker ion intensity}
}
\description{
\code{.PTM_MarkerFinder_no_peakplot} a helper function 
for no peak plot option
}
