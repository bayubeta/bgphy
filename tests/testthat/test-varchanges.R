test_that("prior_transform works", {
  # unbounded stays unbounded
  expect_silent(normal_t <- prior_transform(prior_normal()))
  expect_s3_class(normal_t, c("priorpdf", "normal"))
  expect_equal(attr(normal_t, "bounds"), c(-Inf, Inf))
  expect_equal(attr(normal_t, "param"), stats::setNames(0:1, c("mean", "sd")))

  # lower and upper bounded
  expect_silent(unif_t <- prior_transform(prior_uniform(min = -1, max = 0.5)))
  expect_equal(attr(unif_t, "param"), stats::setNames(c(-1, 0.5), c("min", "max")))
  expect_equal(unif_t(0.2),
               log(dunif(-1 + 1.5*invlogit(0.2), min = -1, max = 0.5)*(1.5)*invlogit(0.2)*(1 - invlogit(0.2))))

  # lower bounded
  expect_silent(gamma_t <- prior_transform(prior_gamma(shape = 2, rate = 3)))
  expect_equal(attr(gamma_t, "param"), stats::setNames(c(2, 3), c("shape", "rate")))
  expect_equal(gamma_t(2.6),
               dgamma(exp(2.6), shape = 2, rate = 3, log = TRUE) + (2.6))
})
