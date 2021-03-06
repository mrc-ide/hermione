% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihood_fast.R
\name{with_recall_bias}
\alias{with_recall_bias}
\title{Probability of observing serial interval taking into account
recall bias}
\usage{
with_recall_bias(tobs, nu, inf_params, ip_params, beta)
}
\arguments{
\item{tobs}{observed serial interval}

\item{nu}{isolation}

\item{inf_params}{parameters of the infectious period}

\item{ip_params}{parameters of the incubation period}

\item{beta}{recall bias coefficient}
}
\value{
probability of observing tobs taking into account recall
bias
}
\description{
Probability of observing serial interval taking into account
recall bias
}
\author{
Sangeeta Bhatia
}
