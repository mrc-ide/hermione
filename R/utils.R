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
