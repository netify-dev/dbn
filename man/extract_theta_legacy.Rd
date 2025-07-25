% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_predict.R
\name{extract_theta_legacy}
\alias{extract_theta_legacy}
\title{Extract theta from legacy format}
\usage{
extract_theta_legacy(fit, draw_idx)
}
\arguments{
\item{fit}{Model fit object}

\item{draw_idx}{Draw index}
}
\value{
Theta array or NULL
}
\description{
Attempt to reconstruct theta array from legacy storage format
}
\keyword{internal}
