% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models_hmm.R
\name{step_ll}
\alias{step_ll}
\title{Step Log-likelihood}
\usage{
step_ll(A, B, Theta_prev, Theta_curr, sigma2_proc)
}
\arguments{
\item{A}{Sender effects matrix}

\item{B}{Receiver effects matrix}

\item{Theta_prev}{Previous Theta}

\item{Theta_curr}{Current Theta}

\item{sigma2_proc}{Innovation variance}
}
\value{
Log-likelihood value
}
\description{
Computes log-likelihood for one time step given regime
}
\keyword{internal}
