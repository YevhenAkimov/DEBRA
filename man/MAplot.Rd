% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{MAplot}
\alias{MAplot}
\title{MA-plot from count means and log2 fold changes}
\usage{
MAplot(drb, FDR = 0.25, filtered = T)
}
\arguments{
\item{drb}{DEBRADataSet object}

\item{FDR}{numeric specifying FDR level. Barcodes above it will be colored in red. Default is 0.25}

\item{filtered}{logical specifying if filtered results will be used to generate MA-plot}
}
\description{
A function that creates MA-plot from DEBRA analysis results
}
\examples{
MAplot(drb, FDR = 0.25)
}
