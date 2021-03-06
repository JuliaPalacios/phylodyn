% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Phyloinfer.R
\name{updateTimes}
\alias{updateTimes}
\title{HMC for times}
\usage{
updateTimes(
  currentval,
  times_list,
  theta_list,
  prec_list,
  probs_list,
  const = 1,
  nsites = nsites,
  eps = 0.04
)
}
\arguments{
\item{times_list}{of coalescent times}

\item{theta_list}{value of log Ne(t)}

\item{prec_list}{vector recording precision parameter}

\item{probs_list}{matrix with likelihod posterior prior at each iteration}

\item{const}{scaling factor}

\item{nsites}{mutation rate}

\item{eps}{leapfrof parameter}

\item{currenval}{current state of the chain}
}
\value{
A list with initial values.
}
\description{
HMC for times
}
