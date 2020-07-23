##' Reparameterise Beta distributions
##'
##'
##'
##' @param mu Mean, should be less than 1
##' @param sigma2 variance, should be less than 1/12
##' @return A named list containing 'shape1' and 'shape2'
##' @author Sangeeta Bhatia
##' @export
beta_muvat2shape1shape2 <- function(mu, sigma2) {

  shape1 <- (mu^2 * (1 - mu) / sigma2) - mu
  shape2 <- shape1 * (1 - mu) / mu
  list(shape1 = shape1, shape2 = shape2)
}

##' Reparameterise Beta distributions
##'
##'
##'
##' @param shape1 positive
##' @param shape2 positive
##' @return A named list containing 'mu' and 'sigma2'
##' @author Sangeeta Bhatia
##' @export
beta_shape1shape22muvar <- function(shape1, shape2) {

  mu <- shape1 / (shape1 + shape2)
  sigma2 <- (shape1 * shape2) / ((shape1 + shape2)^2 * (shape1 + shape2 + 1))
  list(mu = mu, sigma2 = sigma2)
}

##' Modified beta probability density function
##'
##' This function returns the probability density of the random
##' variable x/xmax where x/xmax is distributed according to a beta
##' distribution.
##' @title Modified beta distribution
##' @param x numeric value for which density is required
##' @param shape1 non-negative parameter of the Beta distribution
##' @param shape2 non-negative parameter of the Beta distribution
##' @param xmax Maximum possible value of x
##' @param ncp non-centrality parameter.
##' @param log logical; if TRUE, probabilities p are given as log(p)
##' @return modified density
##' @author Sangeeta Bhatia
##' @export
dmbeta <- function(x, shape1, shape2, xmax, ncp = 0, log = FALSE) {

  dbeta(x = x/xmax, shape1 = shape1, shape2 = shape2, log = log)

}

##' Modififed beta cumulative probability distribution function
##' @inheritParams dmbeta
##' @return modified distribution
##' @author Sangeeta Bhatia
##' @export
pmbeta <- function(x, shape1, shape2, xmax, ncp = 0, log = FALSE) {

  pbeta(x = x/xmax, shape1 = shape1, shape2 = shape2, log = log)

}
