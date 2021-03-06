% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihood.R
\name{probability_isolation_offset}
\alias{probability_isolation_offset}
\title{Probability of serial interval when primary case has been
isolated and in the presence of asymptomatic infectiousness}
\usage{
probability_isolation_offset(
  t,
  nu,
  offset,
  inf_distr = "dgamma",
  inc_distr = "dgamma",
  inf_params,
  inc_params
)
}
\arguments{
\item{t}{observed serial interval}

\item{nu}{onset to isolation of primary case}

\item{offset}{days of asymptomatic infectiousness}

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
numeric probability of observing the given serial interval
with the given parameters of infectious period and incubation
period distributions
}
\description{
Probability of serial interval when primary case has been
isolated and in the presence of asymptomatic infectiousness
}
\author{
Sangeeta Bhatia
}
