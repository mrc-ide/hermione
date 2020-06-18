## Probability of observing a serial interval t is the convolution of
## the distributions of infectious profile and incubation period
##'
##'
##'
##' @title Probability of serial interval as a convolution of
##' infectiousness distribution of primary case and incubation period
##' of secondary case
##' @param t observed serial interval
##' @param tmax the maximum time (from the onset of symptoms in the primary case)
##' at which secondary case could have been infected
##'
##' @param inf_params list with components rate and shape for
##' infectious period distribution
##' @param ip_params list with components rate and shape for
##' incubation period distribution
##' @return numeric. Probability of observing the given serial
##' interval given
##' the parameters of infectious period distributions
##' and incubation period distributions
##' @author Sangeeta Bhatia
##' @export
fast_probability_basic <- function(t, tmax, inf_params, ip_params) {

  f <- function(s) {
    stats::dgamma(s, rate = inf_params$rate, shape = inf_params$shape) *
      stats::dgamma(t - s, rate = ip_params$rate, shape = ip_params$shape)
  }
  out <- stats::integrate(f, 0, tmax)
  out$value
}

##' @author Sangeeta Bhatia
##'
##' @param pbasic Probability without taking into account isolation
##' computed through fast_probability_basic
##' @inheritParams fast_probability_basic
##' @export
fast_probability_isolation <- function(pbasic, nu, inf_params) {

  pbasic /
  stats::pgamma(nu, rate = inf_params$rate, shape = inf_params$shape)
}

## Probability of observing a given serial interval taking into
## account asymptomatic infectiousness.
##' @title Probability of serial interval
##' in the presence of asymptomatic infectiousness
##' @param t observed serial interval
##' @param tmax Max of observed serial interval and isolation if
##' prmary case has been isolated
##' @param offset days of asymptomatic infectiousness
##' @inheritParams probability_basic
##' @return numeric probability of observing the given serial interval
##' with the given parameters of infectious period and incubation
##' period distribution
##' @author Sangeeta Bhatia
##' @export
##'
fast_probability_offset <- function(t, tmax, offset, inf_params, ip_params) {

  f <- function(s) {
    stats::dgamma(
      s + offset,
      rate = inf_params$rate, shape = inf_params$shape
    ) * stats::dgamma(
      t - s,
      rate = ip_params$rate, shape = ip_params$shape
    )
  }

  out <- stats::integrate(f, -offset, tmax)

  out$value
}


## Probability of observing a given serial interval taking into
## account isolation of primary case.
##' @title Probability of serial interval when primary case has been
##' isolated and in the presence of asymptomatic infectiousness
##' @param poffset Probability calculation with offset only
##' @param nu  onset to isolation of primary case
##' @param offset days of asymptomatic infectiousness
##' @inheritParams probability_basic
##' @return numeric probability of observing the given serial interval
##' with the given parameters of infectious period and incubation
##' period distribution
##' @author Sangeeta Bhatia
##' @export
##'
fast_probability_isolation_offset <- function(poffset, nu, offset, inf_params) {

  denominator <- stats::pgamma(
    nu + offset, rate = inf_params$rate, shape = inf_params$shape
  )
  poffset / denominator
}



##' @title Probability of observing serial interval taking into account
##' recall bias
##' @param tobs observed serial interval
##' @param nu isolation
##' @param inf_params  parameters of the infectious period
##' @param ip_params parameters of the incubation period
##' @param beta recall bias coefficient
##' @return probability of observing tobs taking into account recall
##' bias
##' @author Sangeeta Bhatia
##' @export
with_recall_bias <- function(tobs, nu, inf_params, ip_params, beta) {

  f1 <- function(s, tvary) {
    stats::dgamma(s, rate = inf_params$rate, shape = inf_params$shape) *
      stats::dgamma(
        tvary - s, rate = ip_params$rate, shape = ip_params$shape
      )
  }

  f2 <- function(tvary) {
    upper_lim <- min(tvary, nu)
    out <- stats::integrate(f1, lower = 0, upper = upper_lim, t = tvary)
    stats::dexp(tvary, rate = beta) * out$value
  }

  denominator <- stats::pgamma(
    nu, rate = inf_params$rate, shape = inf_params$shape
  ) ##* stats::pexp(tobs, rate = beta)

  out <- stats::integrate(f2, 0, tobs)
  out$value / denominator
}
