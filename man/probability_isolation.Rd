% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihood.R
\name{probability_isolation}
\alias{probability_isolation}
\title{Probability of serial interval}
\usage{
probability_isolation(
  t,
  nu,
  inf_distr = "dgamma",
  inc_distr = "dgamma",
  inf_params,
  inc_params
)
}
\arguments{
\item{t}{observed serial interval}

\item{nu}{onset to isolation of primary case}

\item{inf_distr}{quoted name of the infectious period
density function. Defaults to dgamma}

\item{inc_distr}{quoted name of the incubation period
density function. Defaults to dgamma}

\item{inf_params}{named list of arguments for
infectious period distribution.}

\item{inc_params}{named list of arguments for
incubation period distribution.}
}
\value{
numeric. Probability of observing the given serial
interval given the delay from onset to isolation,
the parameters of infectious period distributions
and incubation period distributions.
}
\description{
Probability of serial interval
}
\author{
Sangeeta Bhatia
}
