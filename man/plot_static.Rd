% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_methods.R
\name{plot_static}
\alias{plot_static}
\title{Plot Static DBN Results}
\usage{
plot_static(results, alpha = 0.01)
}
\arguments{
\item{results}{Output from dbn_static()}

\item{alpha}{Significance level for edge detection (default 0.01)}
}
\value{
A ggplot2 object with multiple panels
}
\description{
Creates diagnostic plots for static model results using ggplot2
}
\keyword{internal}
