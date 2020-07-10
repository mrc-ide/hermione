test_that("probability as convolution of infectious period and incubation period works", {

  alpha1 <- 100
  alpha2 <- 50
  rate <- 100
  p <- probability_basic(
    10,
    inf_params = list(shape = alpha1, rate = rate),
    inc_params = list(shape = alpha2, rate = rate)
  )
  ## convolution of two gamma distributions with the same rate
  ## parameter is a gamma distribution with shape shape1 + shape2
  right <- dgamma(10, shape = 150, rate = 100)
  right <- round(log(p), 4)

  log_p1 <- round(log(p), 4)
  expect_equal(log_p1, right)
})
