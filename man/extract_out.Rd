% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{extract_out}
\alias{extract_out}
\title{Extract posterior beta and effective population size from a stan object and puts them
into the phylodyn format in order to plot}
\usage{
extract_out(
  mfit,
  grid,
  midpts,
  samp_times,
  coal_times,
  n_sampled,
  alpha = 0.05
)
}
\arguments{
\item{mfit}{stan fit object (output of adapref_sampling)}

\item{grid}{grid boundaries}

\item{midpts}{grid mid points}

\item{alpha}{which BCI to keep}
}
\value{
BNPR style output
}
\description{
Extract posterior beta and effective population size from a stan object and puts them
into the phylodyn format in order to plot
}
