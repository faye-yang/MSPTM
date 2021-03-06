% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PTM_marker_finder.R
\name{PTM_MF}
\alias{PTM_MF}
\title{\code{PTM_MF} output graph of .
A helper function of intensity_plot}
\usage{
PTM_MF(data, modification, modi_name, mZmarkerIons, minNumberIons = 2,
  itol_ppm = 10, minMarkerIntensityRatio = 5, mgfFilename = -1,
  PEAKPLOT = TRUE)
}
\arguments{
\item{data}{mass spetrometry information for the peptide}

\item{modification}{contain modification information , intensity of ion, amino acide that is modified}

\item{modi_name}{name of modification}

\item{mZmarkerIons}{maker ion}

\item{minNumberIons}{minimum number or marker ion}

\item{itol_ppm}{ppm}

\item{minMarkerIntensityRatio}{minimum ratio for marker ion intensity}

\item{mgfFilename}{mgf file name indicator}

\item{PEAKPLOT}{peak plot for MTP}
}
\description{
\code{PTM_MF} output graph of .
A helper function of intensity_plot
}
