test_that("Throw an error for invalid input types", {
  expect_error(prior_uniform(min = 0, max = 0))
  expect_error(prior_uniform(min = 0, max = "a"))

  expect_error(prior_normal(mean = 0, sd = -1))
  expect_error(prior_normal(mean = 0, sd = "sd"))

  expect_error(prior_gamma(shape = 0, rate = 1))
  expect_error(prior_gamma(shape = "shape", rate = 1))

  expect_error(prior_halfnormal(sigma = -1))
  expect_error(prior_halfnormal(sigma = "s"))

  expect_error(prior_halfcauchy(sigma = -1))
  expect_error(prior_halfcauchy(sigma = "s"))

  expect_error(prior_halft(nu = 1, sigma = -1))
  expect_error(prior_halft(nu = -1, sigma = 1))
  expect_error(prior_halft(nu = 1, sigma = "s"))
})
