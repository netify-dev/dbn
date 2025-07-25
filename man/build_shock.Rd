% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/irf.R
\name{build_shock}
\alias{build_shock}
\title{Build shock matrix for IRF analysis}
\usage{
build_shock(
  m,
  type = c("unit_edge", "node_out", "node_in", "density"),
  i = 1,
  j = 2,
  magnitude = 1
)
}
\arguments{
\item{m}{Number of nodes in the network}

\item{type}{Type of shock: "unit_edge", "node_out", "node_in", or "density"}

\item{i}{Source node index (for unit_edge and node shocks)}

\item{j}{Target node index (for unit_edge shock)}

\item{magnitude}{Shock magnitude (default: 1)}
}
\value{
m x m shock matrix
}
\description{
Creates shock matrices for different types of network interventions
}
