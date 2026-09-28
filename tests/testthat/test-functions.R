test_that("operator and tensor_vec preserve expected dimensions", {
  X <- diag(3)
  Y <- 1:3

  expect_equal(as.vector(fungarch:::operator(X, Y)), Y / 3)
  expect_equal(dim(fungarch:::tensor_vec(1:3, 4:6)), c(1L, 9L))
  expect_equal(
    as.vector(fungarch:::tensor_vec(1:3, 4:6)),
    as.vector(t(as.vector(1:3) %o% as.vector(4:6)))
  )
})

test_that("operator and tensor_vec reject invalid inputs", {
  expect_error(fungarch:::operator(diag(2), 1:3), "match")
  expect_error(fungarch:::operator(diag(2), c(1, NA)), "finite")
  expect_error(fungarch:::tensor_vec(1:2, c(1, NA)), "finite")
})

test_that("return transformations compute log changes", {
  dat <- matrix(c(1, 2, 2, 4, 4, 8), nrow = 2)

  expect_equal(
    fChange::cidr(dat)$data,
    rbind(rep(0, 3), rep(100 * log(2), 3))
  )
  expect_equal(
    fChange::ocidr(dat)$data,
    matrix(c(0, 100 * log(2), 0, 100 * log(2), 0, 100 * log(2)), nrow = 2)
  )
})

test_that("generateFARCH_use returns simulated data with consistent dimensions", {
  resolution <- 3
  N <- 5
  C_eps <- diag(resolution)
  basis <- diag(sqrt(resolution), resolution)
  alpha <- matrix(0, nrow = resolution^2, ncol = resolution^2)
  delta <- C_eps

  result <- fungarch:::generateFARCH_use(
    N = N,
    resolution = resolution,
    alpha = alpha,
    delta = delta,
    basis = basis,
    burnin = 2,
    p = 1,
    eigs = rep(1, resolution),
    C_eps = C_eps,
    kl_errors = TRUE
  )

  expect_equal(dim(result$data), c(resolution, N))
  expect_equal(dim(result$error_matrix), c(resolution, N))
  expect_equal(length(result$sigma_mat), resolution^2 * N)
})

test_that("generateFARCH_use supports multiple lags", {
  resolution <- 3
  N <- 4
  C_eps <- diag(resolution)
  basis <- diag(sqrt(resolution), resolution)
  alpha_zero <- matrix(0, nrow = resolution^2, ncol = resolution^2)

  result <- fungarch:::generateFARCH_use(
    N = N,
    resolution = resolution,
    alpha = list(alpha_zero, alpha_zero),
    delta = C_eps,
    basis = basis,
    burnin = 2,
    p = 2,
    eigs = rep(1, resolution),
    C_eps = C_eps,
    kl_errors = TRUE
  )

  expect_equal(dim(result$data), c(resolution, N))
})

test_that("new diagonal prediction helper returns square matrices", {
  resolution <- 3
  basis <- diag(sqrt(resolution), resolution)
  C_eps <- diag(resolution)
  delta <- diag(resolution)
  data <- matrix(1, nrow = resolution, ncol = 2)
  alphas <- matrix(c(0.1, 0.2, 0.3), nrow = resolution, ncol = 1)

  prediction <- fungarch:::.pred_diag_new(
    resolution = resolution,
    alphas = alphas,
    basis = basis,
    delta = delta,
    data = data,
    C_eps = C_eps
  )

  expect_equal(dim(prediction), c(resolution, resolution))
  expect_true(all(is.finite(prediction)))
})

test_that("ARCH estimation helper returns parameter estimates", {
  C_0 <- diag(2)
  C_1 <- diag(c(0.2, 0.1))
  projector <- diag(2)
  m2 <- c(1, 2)

  result <- fungarch:::.estimateARCH(
    N = 20,
    n = 2,
    C_0 = C_0,
    C_1 = C_1,
    lambda = 0.1,
    projector = projector,
    m2 = m2
  )

  expect_named(result, c("hat_alpha", "hat_delta"))
  expect_equal(dim(result$hat_alpha), c(2L, 2L))
  expect_equal(length(result$hat_delta), 2L)
})

