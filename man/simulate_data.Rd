% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Tajima.R
\name{simulate_data}
\alias{simulate_data}
\title{Simulate data according to infinite sites mutation model}
\usage{
simulate_data(mu, tree1)
}
\arguments{
\item{mu}{overall mutation rate. The expected number of mutations
will be the (mu)(L) where L is the tree length}

\item{tree1}{ape tree object}
}
\value{
data
}
\description{
Simulate data according to infinite sites mutation model
}
