## Probability of observing a given serial interval taking into
## account isolation of primary case.
##' @title Probability of serial interval when primary case has been
##' isolated
##' @param t observed serial interval
##' @param nu  onset to isolation of primary case
##' @param offset days of asymptomatic infectiousness
##'
##' @return numeric probability of observing the given serial interval
##' with the given parameters of infectious period and incubation
##' period distribution
##' @author Sangeeta Bhatia
##' @export
##'
probability_isolation_offset <- function(t, nu, offset, inf_params, ip_params) {

  f <- function(s) {
    dgamma(
      s + offset,
      rate = inf_params$rate, shape = inf_params$shape
    ) * dgamma(
      t - s,
      rate = ip_params$rate, shape = ip_params$shape
    ) / pgamma(
      nu + offset,
      rate = inf_params$rate, shape = inf_params$shape
    )
  }

  upper_lim <- min(t, nu)
  out <- stats::integrate(f, -offset, upper_lim)

  out$value
}

## Probability of observing a given serial interval taking into
## account isolation of primary case.
##'
##' @param nu onset to isolation of primary case
##' @inheritParams probability_basic
##'
##' @return numeric. Probability of observing the given serial
##' interval given the delay from onset to isolation,
##' the parameters of infectious period distributions
##' and incubation period distributions.
##' @author Sangeeta Bhatia
probability_isolation <- function(t, nu, inf_params, ip_params) {

  f <- function(s) {
    dgamma(s, rate = inf_params$rate, shape = inf_params$shape) *
      dgamma(t - s, rate = ip_params$rate, shape = ip_params$shape) /
      pgamma(nu, rate = inf_params$rate, shape = inf_params$shape)
  }

  upper_lim <- min(t, nu)
  out <- integrate(f, 0, upper_lim)


  out$value
}


## Probability of observing a serial interval t is the convolution of
## the distributions of infectious profile and incubation period
##'
##'
##'
##'
##' @param t observed serial interval
##' @param inf_params
##' @param ip_params
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
probability_basic <- function(t, inf_params, ip_params) {

  f <- function(s) {
    dgamma(s, rate = inf_params$rate, shape = inf_params$shape) *
      dgamma(t - s, rate = ip_params$rate, shape = ip_params$shape)
  }
  out <- stats::integrate(f, 0, t)
  out$value
}

## Provide a named argument to this function
## For example:
## log_likelihood(t, inf_params, ip_params, probability_isolation, nu = x)
log_likelihood <- function(t, inf_params, ip_params, fun, ...) {
  log(fun(t, inf_params, ip_params, ...))
}

total_log_likelihood <- function(tvec,
                                 nu_vec = NULL,
                                 offset_vec = NULL,
                                 inf_params,
                                 ip_params) {
  if ((!is.null(offset_vec)) & (!is.null(nu_vec))) {
    ## Use probability_isolation_offset
    out <- mapply(
      FUN = log_likelihood,
      t = tvec,
      nu = nu_vec,
      offset = offset_vec,
      MoreArgs = list(
        inf_params = inf_params,
        ip_params = ip_params,
        fun = probability_isolation_offset
      ),
      SIMPLIFY = TRUE
    )
  } else if (!is.null(nu_vec)) {
    ## Use probability_isolation
    out <- mapply(
      FUN = log_likelihood,
      t = tvec,
      nu = nu_vec,
      MoreArgs = list(
        inf_params = inf_params,
        ip_params = ip_params,
        fun = probability_isolation
      ),
      SIMPLIFY = TRUE
    )
  } else {
    ## Use probability_basic
    out <- sapply(
      tvec,
      function(t) log_likelihood(t, inf_params, ip_params, fun = probability_basic)
    )
  }

  # out <- out[is.finite(out)]
  sum(out)
}
