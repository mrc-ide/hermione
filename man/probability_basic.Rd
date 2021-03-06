% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihood.R
\name{probability_basic}
\alias{probability_basic}
\title{Probability of serial interval as a convolution of
infectiousness distribution of primary case and incubation period
of secondary case}
\usage{
probability_basic(
  t,
  inf_distr = "dgamma",
  inc_distr = "dgamma",
  inf_params,
  inc_params
)
}
\arguments{
\item{t}{observed serial interval}

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
interval given
the parameters of infectious period distributions
and incubation period distributions
}
\description{
Probability of serial interval
}
\details{
Probability of observing a serial interval t is the convolution of
the distributions of infectious profile and incubation period.
}
\examples{
probability_basic(
  10,
  inf_params = list(shape = 100, rate = 100),
  inc_params = list(shape = 50, rate = 100)
)
}
\author{
Sangeeta Bhatia
}
