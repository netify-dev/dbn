% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/irf.R
\name{compute_irf_single}
\alias{compute_irf_single}
\title{Compute IRF for a single posterior draw}
\usage{
compute_irf_single(fit, draw_idx, shock, H, t0 = 1, stat_fun = stat_density)
}
\arguments{
\item{fit}{dbn model fit object}

\item{draw_idx}{Index of posterior draw to use}

\item{shock}{Shock matrix}

\item{H}{Number of horizons}

\item{t0}{Shock time (for dynamic models, 1-based)}

\item{stat_fun}{Network statistic function}
}
\value{
Vector of IRF values at each horizon
}
\description{
Compute IRF for a single posterior draw
}
\keyword{internal}
