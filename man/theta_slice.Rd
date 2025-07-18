% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_extract.R
\name{theta_slice}
\alias{theta_slice}
\title{Extract Theta slices from posterior draws}
\usage{
theta_slice(fit, draws = NULL, i = NULL, j = NULL, rel = NULL, time = NULL)
}
\arguments{
\item{fit}{A dbn model fit object}

\item{draws}{Integer vector of draw indices to extract}

\item{i}{Row indices (sender nodes)}

\item{j}{Column indices (receiver nodes)}

\item{rel}{Relation indices}

\item{time}{Time indices}
}
\value{
List of Theta slices
}
\description{
Extract specific slices of Theta arrays from posterior draws
}
\keyword{internal}
