##' Simulate serial interval
##'
##' This function simulates serial interval with or without case
##' isolation, and with or without pre-symptomatic infectiousness.
##' At the moment, infectious period, incubation period and onset to
##' isolation distributions are assumed to be gamma distributed.
##' @title Simulate serial interval
##' @param mean_ip Mean Incubation Period
##' @param sd_ip  Standard deviation of the incubation period
##' @param mean_inf Infectious Period
##' @param sd_inf Standard deviation of the infectious period
##' @param mean_iso Mean delay from symptom onset to case isolation.
##' If this is null, then isolation is not simulated.
##' @param sd_iso Standard deviation of the delay from symptom onset
##' to case isolation
##'
##' @param offset Numeric - days before symptom onset that
##' infectiousness is assumed to begin. If \code{offset} is NULL, then
##' asymptomatic infectiousness is assumed to be absent.
##' @param nsim Number of transmission pairs to simulate
##' @return a data frame with columns t_1, t_2, si
##' t_1 corresponds to draws from the infectious period,
##' t_2 corresponds to draws from the incubation period, and
##' si is the serial interval (t_1 + t_2).
##' If mean_iso and sd_iso are not NULL, then the output also contains
##' a column nu corresponding to draws from the distribution of onset
##' to isolation.
##' @author Sangeeta Bhatia
##' @importFrom epitrix gamma_mucv2shapescale
##' @export
##'
simulate_si <- function(mean_ip,
                        sd_ip,
                        mean_inf,
                        sd_inf,
                        mean_iso = NULL,
                        sd_iso = NULL,
                        offset = 0,
                        nsim = 50) {
  params_ip <- epitrix::gamma_mucv2shapescale(
    mu = mean_ip, cv = sd_ip / mean_ip
  )

  mean_inf <- mean_inf + offset

  params_inf <- epitrix::gamma_mucv2shapescale(
    mu = mean_inf, cv = sd_inf / mean_inf
  )

  t_1 <- stats::rgamma(
    n = nsim, shape = params_inf$shape, rate = 1 / params_inf$scale
  )

  ## We shifted the mean of the distribution to the right
  ## And then shift back the samples by amount offset.
  t_1 <- t_1 - offset
  t_2 <- stats::rgamma(
    n = nsim, shape = params_ip$shape, rate = 1 / params_ip$scale
  )

  out <- data.frame(t_1 = t_1, t_2 = t_2, si = round(t_1 + t_2))

  if (!is.null(mean_iso)) {

    params_iso <- epitrix::gamma_mucv2shapescale(
      mu = mean_iso, cv = sd_iso / mean_iso
    )
    out$nu <- stats::rgamma(
      n = nsim, shape = params_iso$shape, rate = 1 / params_iso$scale
    )
    out$nu <- round(out$nu)
    ## t_1 i.e. the time at which secondary could have been infected
    ## must be less than nu which is the time from symptom onset to
    ## isolation.
    out <- out[out$t_1 < out$nu, ]
    ## TODO Make sure you get the desired number of rows after
    ## applying these filters.
  }
  ## 0 Serial Interval is not allowed if we do not account for
  ## pre-symptomatic infectivity.
  if (offset == 0) out <- out[out$si > 0, ]

  out
}
