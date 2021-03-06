## Probability of observing a given serial interval taking into
## account asymptomatic infectiousness.
##' @title Probability of serial interval
##' with asymptomatic infectiousness
##' @param t observed serial interval
##'
##' @param offset days of asymptomatic infectiousness
##'
##' @inheritParams probability_basic
##' @return numeric probability of observing the given serial interval
##' with the given parameters of infectious period and incubation
##' period distribution
##' @author Sangeeta Bhatia
##' @export
##'
probability_offset <- function(t,
                               offset,
                               inf_distr = "dgamma",
                               inc_distr = "dgamma",
                               inf_params,
                               inc_params) {

  f1 <- match.fun(inf_distr)
  f2 <- match.fun(inc_distr)
  f <- function(s) {
    inf_params$x <- s + offset
    inc_params$x <- t - s
    do.call(f1, inf_params) * do.call(f2, inc_params)
  }

  out <- stats::integrate(f, -offset, t, stop.on.error = FALSE)

  out$value
}


## Probability of observing a given serial interval taking into
## account isolation of primary case.
##' @title Probability of serial interval when primary case has been
##' isolated and in the presence of asymptomatic infectiousness
##' @param t observed serial interval
##' @param nu onset to isolation of primary case
##' @param offset days of asymptomatic infectiousness
##'
##'
##'
##'
##' @inheritParams probability_basic
##' @return numeric probability of observing the given serial interval
##' with the given parameters of infectious period and incubation
##' period distributions
##' @author Sangeeta Bhatia
##' @export
##'
probability_isolation_offset <- function(t,
                                         nu,
                                         offset,
                                         inf_distr = "dgamma",
                                         inc_distr = "dgamma",
                                         inf_params,
                                         inc_params) {


  f1 <- match.fun(inf_distr)
  f2 <- match.fun(inc_distr)
  f <- function(s) {
    inf_params$x <- s + offset
    inc_params$x <- t - s
    do.call(f1, inf_params) * do.call(f2, inc_params)
  }

  upper_lim <- min(t, nu)
  out <- stats::integrate(f, -offset, upper_lim, stop.on.error = FALSE)

  g <- gsub("^d", "p", inf_distr)
  inf_params$q <- nu
  denominator <- do.call(
    what = g,
    args = inf_params
  )
  out$value / denominator
}

## Probability of observing a given serial interval taking into
## account isolation of primary case.
##' @title Probability of serial interval
##'
##' @param nu onset to isolation of primary case
##' @inheritParams probability_basic
##'
##' @return numeric. Probability of observing the given serial
##' interval given the delay from onset to isolation,
##' the parameters of infectious period distributions
##' and incubation period distributions.
##' @author Sangeeta Bhatia
##' @export
probability_isolation <- function(t,
                                  nu,
                                  inf_distr = "dgamma",
                                  inc_distr = "dgamma",
                                  inf_params,
                                  inc_params) {

  f1 <- match.fun(inf_distr)
  f2 <- match.fun(inc_distr)
  f <- function(s) {
    inf_params$x <- s
    inc_params$x <- t - s
    do.call(f1, inf_params) * do.call(f2, inc_params)
  }

  upper_lim <- min(t, nu)
  out <- stats::integrate(f, 0, upper_lim, stop.on.error = FALSE)
  ## Get the name of fumulative density function
  ## from the name of the density function
  g <- gsub("^d", "p", inf_distr)
  inf_params$q <- nu
  denominator <- do.call(what = g, args = inf_params)

  out$value / denominator
}


##' Probability of serial interval
##'
##' @details
##' Probability of observing a serial interval t is the convolution of
##' the distributions of infectious profile and incubation period.
##'
##'
##'
##' @title Probability of serial interval as a convolution of
##' infectiousness distribution of primary case and incubation period
##' of secondary case
##' @param t observed serial interval
##'
##' @param inf_distr quoted name of the infectious period
##' density function. Defaults to dgamma
##' @param inc_distr quoted name of the incubation period
##' density function. Defaults to dgamma
##' @param inf_params named list of arguments for
##' infectious period distribution.
##' @param inc_params named list of arguments for
##' incubation period distribution.
##'
##'
##' @return numeric. Probability of observing the given serial
##' interval given
##' the parameters of infectious period distributions
##' and incubation period distributions
##' @author Sangeeta Bhatia
##' @export
##' @examples
##' probability_basic(
##'   10,
##'   inf_params = list(shape = 100, rate = 100),
##'   inc_params = list(shape = 50, rate = 100)
##' )
probability_basic <- function(t,
                              inf_distr = "dgamma",
                              inc_distr = "dgamma",
                              inf_params,
                              inc_params) {

  f1 <- match.fun(inf_distr)
  f2 <- match.fun(inc_distr)
  f <- function(s) {
    inf_params$x <- s
    inc_params$x <- t - s
    do.call(f1, inf_params) * do.call(f2, inc_params)
  }
  out <- stats::integrate(f, 0, t, stop.on.error = FALSE)
  out$value
}

## Provide a named argument to this function
## For example:
## log_likelihood(t, inf_params, ip_params, probability_isolation, nu = x)
log_likelihood <- function(t, inf_params, ip_params, fun, ...) {
  log(fun(t, inf_params, ip_params, ...))
}

##' Total log-likelhood
##'
##' @details
##' This function returns the totoal log-likelihood of the  observed
##' serial interval. In the absence of isolation, the likelihood is
##' \deqn{L(rate, shape \mid tvec) = \prod\limits_{i = 1}{n}{ = \int\limits_{0}^{t_{i}}{f(t_1)g(t_{i} - t_1)dt_1} }}
##' where \eqn{t_i} is the ith observation, \eqn{t_1} is the delay between
##' symptom onset in a primary case to the infection of a secondary case,
##' and f and g are the infectious and incubation period distributions
##' respectively.
##' If isolation is included via non-null \code{nu_vec}, the likelihood
##' is
##' \deqn{L(rate, shape \mid tvec) = \prod\limits_{i = 1}{n}{ = \int\limits_{0}^{t_{i}}{(f(t_1) / F(\nu_i))g(t_{i} - t_1)dt_1} }}
##' where \eqn{\nu_i} is the delay from symptom onset in the primary
##' case to their isolation.
##' Finally, if asymptomatic infectiosuness is included (non-null \code{offset_vec})
##' the infectiousness profile is a gamma distribution with offset.
##' @title Total log likelihood
##' @param tvec Vector of observed serial intervals
##' @param nu_vec Vector of delays from onset to isolation
##' @param offset_vec Vector specifying the pre-symptomatic infectiousness
##' @inheritParams probability_basic
##'
##' @return total log-likelihood
##' @author Sangeeta Bhatia
##' @export
total_log_likelihood <- function(tvec,
                                 nu_vec = NULL,
                                 offset_vec = NULL,
                                 inf_params,
                                 ip_params) {

  if ((!is.null(offset_vec)) & (!is.null(nu_vec))) {

    message("Both offset and isolation are non-null")
    message("Using probability_isolation_offset")
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
  } else if (!is.null(nu_vec) & (is.null(offset_vec))) {
    ## Use probability_isolation
    message("Isolation is non-null, offset is NULL")
    message("Using probability_isolation")

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
  } else if (is.null(nu_vec) & (!is.null(offset_vec))) {
    ## Use probability_offset
    message("Isolation is null, offset is not NULL")
    message("Using probability_offset")

    out <- mapply(
      FUN = log_likelihood,
      t = tvec,
      offset = offset_vec,
      MoreArgs = list(
        inf_params = inf_params,
        ip_params = ip_params,
        fun = probability_offset
      ),
      SIMPLIFY = TRUE
    )
  } else if (is.null(nu_vec) & is.null(offset_vec)){
    ## Use probability_basic
    message("Isolation is null, offset is NULL")
    message("Using probability_basic")

    out <- sapply(
      tvec,
      function(t) log_likelihood(t, inf_params, ip_params, fun = probability_basic)
    )
  }

  # out <- out[is.finite(out)]
  sum(out)
}
