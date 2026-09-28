test_that("check_onestep completes on a small simulated dataset", {
  skip("Legacy check_onestep() is currently disabled")
  skip_if_not_installed("fChange")
  skip_on_cran()

  set.seed(5)
  data <- matrix(rnorm(3 * 8), nrow = 3)
  result <- fungarch:::check_onestep(
    data,
    start = 4,
    est_dimen = 2,
    silent = TRUE
  )

  expect_equal(dim(result), c(4L, 7L))
  expect_true(all(is.finite(result)))
})

test_that("compare_models completes on a small simulated dataset", {
  skip("Legacy compare_models() is currently disabled")
  skip_if_not_installed("fChange")
  skip_on_cran()

  set.seed(6)
  result <- fungarch:::compare_models(
    N = 5,
    resolution = 3,
    alpha_scale = c(0.1, 0.1),
    est_dimen = 2,
    N_fit = 3,
    silent = TRUE
  )

  expect_equal(dim(result), c(2L, 9L))
  expect_true(all(is.finite(result)))
})

