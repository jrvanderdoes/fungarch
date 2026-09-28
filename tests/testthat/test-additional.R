test_that("generate_data returns the requested number of observations", {
  set.seed(3)

  result <- fungarch:::generate_data(
    N = 4,
    resolution = 3,
    alpha_scale = 0.1
  )

  expect_equal(dim(result$data), c(3L, 4L))
  expect_equal(dim(result$error_matrix), c(3L, 4L))
  expect_equal(length(result$sigma_mat), 3L^2 * 4L)
})

test_that("pointwise fitting returns coefficient dimensions", {
  set.seed(4)
  X <- matrix(rnorm(3 * 40), nrow = 3)

  result <- fungarch:::.fit_pw(X, inv_type = "M")

  expect_named(result, c("alpha", "delta"))
  expect_equal(dim(result$alpha), c(3L, 3L))
  expect_equal(length(result$delta), 3L)
})

test_that("pointwise fitting rejects invalid inversion types", {
  X <- matrix(rnorm(3 * 10), nrow = 3)

  expect_error(
    fungarch:::.fit_pw(X, inv_type = "X"),
    "one of"
  )
})

test_that("pointwise fitting validates inputs and training splits", {
  X <- matrix(rnorm(3 * 12), nrow = 3)

  expect_error(fungarch:::.fit_pw(X[, 1:2, drop = FALSE]), "at least three")
  expect_error(fungarch:::.fit_pw(X, w = -0.1), "between 0 and 1")
  expect_error(fungarch:::.fit_pw(X, inv_type = "T", training = 0.1), "at least three")
  expect_error(fungarch:::.fit_pw(X, training = 0), "positive finite")
})

test_that("pointwise Tikhonov fitting supports a training split", {
  skip_if_not_installed("fChange")

  set.seed(5)
  X <- matrix(rnorm(3 * 12), nrow = 3)
  result <- fungarch:::.fit_pw(X, inv_type = "T", training = 0.75)

  expect_named(result, c("alpha", "delta"))
  expect_equal(dim(result$alpha), c(3L, 3L))
  expect_true(all(is.finite(result$alpha)))
  expect_true(all(is.finite(result$delta)))
})

test_that("new pointwise fitting returns one operator per lag", {
  set.seed(6)
  X <- matrix(rnorm(3 * 20), nrow = 3)
  result <- fungarch:::.fit_pw_new(X, p = 2)

  expect_length(result$alpha, 2L)
  expect_true(all(vapply(result$alpha, function(x) {
    is.matrix(x) && all(dim(x) == c(3L, 3L)) && all(is.finite(x))
  }, logical(1))))
  expect_length(result$delta, 3L)
  expect_length(result$mean, 3L)
})

test_that("new pointwise prediction uses all requested lags", {
  skip_if_not_installed("fChange")

  X <- matrix(c(1, 2, 2, 4, 3, 6), nrow = 2)
  alpha <- list(diag(c(0.1, 0.2)), diag(c(0.3, 0.4)))
  delta <- c(0.2, 0.3)
  C_eps <- diag(2)

  prediction <- fungarch:::.pred_pw_new(
    X, alpha, delta, C_eps, mean = c(0, 0)
  )
  simple_prediction <- fungarch:::.pred_pw_new_simple(
    X, alpha, delta, mean = c(0, 0)
  )

  expect_equal(dim(prediction), c(2L, 2L))
  expect_length(simple_prediction, 2L)
  expect_true(all(is.finite(prediction)))
  expect_true(all(is.finite(simple_prediction)))
})

test_that("generateFARCH_use validates innovation inputs", {
  resolution <- 3
  basis <- diag(sqrt(resolution), resolution)
  alpha <- matrix(0, nrow = resolution^2, ncol = resolution^2)
  delta <- diag(resolution)

  expect_error(
    fungarch:::generateFARCH_use(
      N = 2,
      resolution = resolution,
      alpha = alpha,
      delta = delta,
      basis = basis,
      kl_errors = TRUE
    ),
    "eigs must be supplied"
  )
  expect_error(
    fungarch:::generateFARCH_use(
      N = 2,
      resolution = resolution,
      alpha = alpha,
      delta = delta,
      basis = basis,
      kl_errors = FALSE
    ),
    "C_eps must be supplied"
  )
})

test_that("data generation validates dimensions and simulation settings", {
  resolution <- 3
  basis <- diag(sqrt(resolution), resolution)
  alpha <- matrix(0, nrow = resolution^2, ncol = resolution^2)
  delta <- diag(resolution)

  expect_error(
    fungarch:::generate_data(N = 0, resolution = resolution),
    "positive integer"
  )
  expect_error(
    fungarch:::generate_data(N = 2, resolution = resolution, alpha_scale = rep(0.1, 4)),
    "between one and resolution"
  )
  expect_error(
    fungarch:::generateFARCH_use(
      N = 2, resolution = resolution, alpha = alpha, delta = delta,
      basis = basis, burnin = 1, p = 4, eigs = rep(1, resolution),
      kl_errors = TRUE
    ),
    "p must be smaller than burnin [+] N"
  )
  expect_error(
    fungarch:::generateFARCH_use(
      N = 2, resolution = resolution, alpha = alpha, delta = delta,
      basis = basis, burnin = 1, p = 1, eigs = rep(-1, resolution),
      kl_errors = TRUE
    ),
    "positive finite"
  )
})

test_that("new diagonal prediction helpers return finite matrices", {
  resolution <- 3
  basis <- diag(sqrt(resolution), resolution)
  C_eps <- diag(resolution)
  delta <- diag(resolution)
  data <- matrix(1, nrow = resolution, ncol = 4)
  alphas <- matrix(0.1, nrow = 2, ncol = 1)

  prediction <- fungarch:::.pred_diag_new(
    resolution = resolution,
    alphas = alphas,
    basis = basis,
    delta = delta,
    data = data,
    C_eps = C_eps
  )
  sigma_sqrt <- fungarch:::.sigma2_diag(
    resolution = resolution,
    alphas = alphas,
    basis = basis,
    delta = delta,
    data = data,
    C_eps = C_eps
  )

  expect_equal(dim(prediction), c(resolution, resolution))
  expect_equal(dim(sigma_sqrt), c(resolution, resolution))
  expect_true(all(is.finite(prediction)))
  expect_true(all(is.finite(sigma_sqrt)))
})

test_that("new diagonal prediction uses lag-specific coefficients", {
  resolution <- 2
  basis <- diag(sqrt(resolution), resolution)
  C_eps <- diag(resolution)
  delta <- matrix(0, nrow = resolution, ncol = resolution)
  data <- matrix(c(1, 2, 2, 4, 3, 6), nrow = resolution)
  alphas <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2, ncol = 2)
  swapped_alphas <- alphas[, c(2, 1)]

  prediction <- fungarch:::.pred_diag_new(
    resolution = resolution,
    alphas = alphas,
    basis = basis,
    delta = delta,
    data = data,
    C_eps = C_eps
  )
  swapped_prediction <- fungarch:::.pred_diag_new(
    resolution = resolution,
    alphas = swapped_alphas,
    basis = basis,
    delta = delta,
    data = data,
    C_eps = C_eps
  )

  expect_true(all(is.finite(prediction)))
  expect_false(isTRUE(all.equal(prediction, swapped_prediction)))
})

test_that("new diagonal fitting falls back to the full basis", {
  skip_if_not_installed("fChange")

  resolution <- 3
  set.seed(7)
  X <- matrix(rnorm(resolution * 20), nrow = resolution)
  C_eps <- diag(resolution)
  basis <- diag(sqrt(resolution), resolution)
  eigenvalues <- rep(1, resolution)

  result <- fungarch:::.fit_diag_model_new(
    X = X,
    C_eps = C_eps,
    basis = basis,
    eigenvalues = eigenvalues,
    w = 0,
    p = 1
  )

  expect_equal(result$K, resolution)
})

test_that("new diagonal fitting validates its inputs", {
  resolution <- 3
  X <- matrix(1, nrow = resolution, ncol = 10)
  C_eps <- diag(resolution)
  basis <- diag(sqrt(resolution), resolution)
  eigenvalues <- rep(1, resolution)

  expect_error(
    fungarch:::.fit_diag_model_new(
      X = X[-1, , drop = FALSE], C_eps = C_eps, basis = basis,
      eigenvalues = eigenvalues
    ),
    "matching X"
  )
  expect_error(
    fungarch:::.fit_diag_model_new(
      X = X, C_eps = C_eps, basis = basis[-1, , drop = FALSE],
      eigenvalues = eigenvalues
    ),
    "matching rows"
  )
  expect_error(
    fungarch:::.fit_diag_model_new(
      X = X, C_eps = C_eps, basis = basis, eigenvalues = eigenvalues,
      p = 0
    ),
    "positive integer"
  )
  expect_error(
    fungarch:::.fit_diag_model_new(
      X = X, C_eps = C_eps, basis = basis, eigenvalues = eigenvalues,
      w = 1.1
    ),
    "between 0 and 1"
  )
})

test_that("new diagonal fitting supports multiple lags", {
  skip_if_not_installed("fChange")

  resolution <- 3
  set.seed(8)
  X <- matrix(rnorm(resolution * 30), nrow = resolution)
  C_eps <- diag(resolution)
  basis <- diag(sqrt(resolution), resolution)
  eigenvalues <- rep(1, resolution)

  result <- fungarch:::.fit_diag_model_new(
    X = X,
    C_eps = C_eps,
    basis = basis,
    eigenvalues = eigenvalues,
    w = 0,
    p = 2
  )

  expect_equal(dim(result$alpha), c(resolution, 2L))
  expect_true(all(is.finite(result$alpha)))
  expect_true(all(is.finite(result$delta)))
})

test_that("new diagonal fitting accepts a smaller supplied basis", {
  skip_if_not_installed("fChange")

  resolution <- 3
  set.seed(9)
  X <- matrix(rnorm(resolution * 20), nrow = resolution)
  C_eps <- diag(resolution)
  basis <- diag(sqrt(resolution), resolution)[, 1:2, drop = FALSE]
  eigenvalues <- rep(1, ncol(basis))

  result <- fungarch:::.fit_diag_model_new(
    X = X,
    C_eps = C_eps,
    basis = basis,
    eigenvalues = eigenvalues,
    w = 0,
    p = 1
  )

  expect_equal(result$K, 2L)
})

test_that("Tikhonov diagonal fitting supports a training split", {
  skip_if_not_installed("fChange")

  resolution <- 2
  set.seed(10)
  X <- matrix(rnorm(resolution * 12), nrow = resolution)
  C_eps <- diag(resolution)
  basis <- diag(sqrt(resolution), resolution)
  eigenvalues <- rep(1, resolution)

  result <- fungarch:::.fit_diag_model_tik(
    X = X,
    C_eps = C_eps,
    basis = basis,
    eigenvalues = eigenvalues,
    w = 0,
    p = 1,
    lambda = 0.5,
    training = 0.75,
    alp = 0.5
  )

  expect_named(result, c("alpha", "delta", "K_Ce"))
  expect_true(all(is.finite(result$alpha)))
  expect_true(all(is.finite(result$delta)))
})

test_that("pointwise prediction helpers return expected shapes", {
  skip_if_not_installed("fChange")

  X <- matrix(1:12, nrow = 3)
  alpha <- diag(0.1, 3)
  delta <- rep(0.2, 3)
  C_eps <- diag(3)

  prediction <- fungarch:::.pred_pw(X, alpha, delta, C_eps)
  simple_prediction <- fungarch:::.pred_pw_simple(X, alpha, delta, C_eps)

  expect_equal(dim(prediction), c(3L, 3L))
  expect_length(simple_prediction, 3L)
  expect_true(all(is.finite(prediction)))
  expect_true(all(is.finite(simple_prediction)))
})

test_that("ARCH prediction and error helpers return finite values", {
  skip_if_not_installed("fChange")

  X2 <- matrix(c(1, 2, 2, 3, 3, 4), nrow = 2)
  alpha <- diag(0.1, 2)
  delta <- c(0.2, 0.2)

  estimates <- fungarch:::.predictARCH(X2, alpha, delta)
  error <- fungarch:::.computeARCHError(X2, alpha, delta)

  expect_equal(dim(estimates), dim(X2))
  expect_length(error, 1L)
  expect_true(all(is.finite(estimates)))
  expect_true(is.finite(error))
})

