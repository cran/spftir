% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spftir.R
\name{spoffs}
\alias{spoffs}
\title{Offset Correction of a Matrix of Mid-infrared Spectra}
\usage{
spoffs(spectra)
}
\arguments{
\item{spectra}{matrix. The matrix of FTIR spectra. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.}
}
\value{
A matrix spectra with with background values deleted. The first row corresponds to wavenumber; the remaining rows corresponds to absorbances.
}
\description{
Allows the removal of the background of a matrix of mid-infrared spectra.
}
\examples{
data(spData)
# Convert data frame to matrix
spectra <- as.matrix(t(spData))
# Offset correction
offs <- spoffs(spectra=spectra)
}
\author{
Claudio Pozo Valenzuela [aut, cre] and Saddys Rodriguez-llamazares [aut]
}

