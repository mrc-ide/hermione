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

  ## convolution of epxonential distributions is gamma
  p <- probability_basic(
    10, "dexp", "dexp", list(rate = 5), list(rate = 5)
  )
  right <- dgamma(10, 2, 5)
  right <- round(log(p), 4)
  log_p1 <- round(log(p), 4)
  expect_equal(log_p1, right)

})


test_that("probability with pre-symptomatic infectiousness works", {

  p <- probability_offset(
    2, 2, "dexp", "dexp", list(rate = 1), list(rate = 1)
  )
  ## Working out the math leads to:
  ## h(t) = integral from -t to offset with respect to s of
  ## exp(-(t + offset)) which is equal to
  ## (t + offset) * exp(-(t + offset)) which is the exponetial
  ## distribution expression with rate = t + offset and x = 1
  right <- dexp(1, rate = 4)
  right <- round(log(p), 4)

  log_p <- round(log(p), 4)
  expect_equal(log_p, right)

})
