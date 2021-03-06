% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{simulateUnivariate}
\alias{simulateUnivariate}
\title{Simulate univariate data}
\usage{
simulateUnivariate(bins, transition, emission, fragLen = 50)
}
\arguments{
\item{bins}{A \code{\link[GenomicRanges]{GRanges}} object for which reads will be simulated.}

\item{transition}{A matrix with transition probabilities.}

\item{emission}{A data.frame with emission distributions (see \code{\link{uniHMM}} entry 'distributions').}

\item{fragLen}{Length of the simulated read fragments.}
}
\value{
A \code{list} with entries $bins containing the simulated states and read count, $reads with simulated read coordinates and $transition and $emission.
}
\description{
Simulate known states, read counts and read coordinates using a univariate Hidden Markov Model with three states ("zero-inflation", "unmodified" and "modified").
}

