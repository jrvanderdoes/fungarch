test_that("white-noise alpha operators are lag-specific", {
  basis <- diag(sqrt(2), 2)
  alphas <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)

  operators <- fungarch:::.build_white_noise_alpha_operators(
    alphas = alphas,
    basis = basis,
    resolution = 2
  )

  h1 <- as.vector(basis[, 1] %o% basis[, 1])
  h2 <- as.vector(basis[, 2] %o% basis[, 2])
  expected_lag_1 <- 1 * (h1 %o% h1) + 3 * (h2 %o% h2)
  expected_lag_2 <- 2 * (h1 %o% h1) + 4 * (h2 %o% h2)

  expect_equal(operators[[1]], expected_lag_1)
  expect_equal(operators[[2]], expected_lag_2)
  expect_false(isTRUE(all.equal(operators[[1]], operators[[2]])))
})
