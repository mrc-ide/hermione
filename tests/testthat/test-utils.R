test_that("Mean and variance to shape1 and shape2 works", {

  mu <- 0.3
  sigma2 <- 0.01
  out <- beta_muvat2shape1shape2(mu, sigma2)
  right <- beta_shape1shape22muvar(out$shape1, out$shape2)
  expect_equal(right$mu, mu)
  expect_equal(right$sigma2, sigma2)
})


test_that("shape1 and shape2 to mean and variance works", {

  shape1 <- 10
  shape2 <- 10
  out <- beta_shape1shape22muvar(shape1, shape2)
  right <- beta_muvat2shape1shape2(out$mu, out$sigma2)
  expect_equal(right$shape1, shape1)
  expect_equal(right$shape2, shape2)
})
