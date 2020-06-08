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
##' @export
fast_probability_isolation <- function(pbasic, nu, inf_params) {
  
  pbasic / 
  stats::pgamma(nu, rate = inf_params$rate, shape = inf_params$shape)
}