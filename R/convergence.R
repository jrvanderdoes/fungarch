#' Compare convergence of functional ARCH parameter estimates
#'
#' Runs repeated simulations and compares the convergence of diagonal ARCH(1)
#' and ARCH(5) parameter estimates. The target can be either the alpha
#' operators or the innovation covariance contribution.
#'
#' @param parameter Target parameter to evaluate: `"alpha"` or `"delta"`.
#' @param Ns Integer vector of sample sizes.
#' @param resolution Number of grid points in each simulated curve.
#' @param nSims Number of simulations for each sample size.
#' @param alpha_scale Weights used to construct the true alpha operator.
#' @param ap_amt Lag-specific alpha magnitudes.
#' @param p True process lag order.
#' @param w Smoothing/regularization weight used in model fitting.
#' @param burnin Number of simulated burn-in observations.
#' @param silent If `TRUE`, suppress progress output.
#' @param parallel If `TRUE`, run simulations with `future.apply`.
#' @param ncores Number of worker processes when `parallel = TRUE`.
#' @param checkpoint_dir Optional directory for per-sample-size checkpoint files.
#' @return A list containing mean `L2`, `Abs`, and `Rel` error matrices, plus
#'   simulation-level error lists.
#' @export
compare_convergence <- function(
    parameter = c("alpha", "delta"),
    Ns = c(50, 100, 250, 500, 750, 1000),
    resolution = 40,
    nSims = 500,
    alpha_scale = 0.7,
    ap_amt = 0.3,
    p = 1,
    w = 0.1,
    burnin = 100,
    silent = FALSE,
    parallel = FALSE,
    ncores = max(1L, parallel::detectCores(logical = TRUE) - 1L),
    checkpoint_dir = NULL) {
  parameter <- match.arg(parameter)

  if (!is.numeric(Ns) || any(Ns < 1) || any(Ns != as.integer(Ns))) {
    stop("Ns must contain positive integers.")
  }
  if (!is.numeric(resolution) || length(resolution) != 1 ||
      resolution != as.integer(resolution) || resolution < 2) {
    stop("resolution must be one integer greater than 1.")
  }
  if (!is.numeric(nSims) || length(nSims) != 1 ||
      nSims != as.integer(nSims) || nSims < 1) {
    stop("nSims must be one positive integer.")
  }
  if (!is.numeric(p) || length(p) != 1 || p != as.integer(p) || p < 1) {
    stop("p must be one positive integer.")
  }
  if (length(ap_amt) == 1 && p > 1) ap_amt <- rep(ap_amt, p)
  if (length(ap_amt) != p) {
    stop("ap_amt must have length one or length p.")
  }
  if (length(alpha_scale) > resolution) {
    stop("alpha_scale cannot contain more basis weights than resolution.")
  }
  if (!is.logical(parallel) || length(parallel) != 1) {
    stop("parallel must be TRUE or FALSE.")
  }
  if (!is.numeric(ncores) || length(ncores) != 1 ||
      ncores != as.integer(ncores) || ncores < 1) {
    stop("ncores must be one positive integer.")
  }
  if (!is.null(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop("parallel = TRUE requires the future and future.apply packages.")
    }
    old_plan <- future::plan()
    future::plan(future::multisession, workers = as.integer(ncores))
    on.exit(future::plan(old_plan), add = TRUE)
  }

  C_eps <- .generate_ceps(resolution)
  eig_Ce <- eigen(C_eps)
  basis <- eig_Ce$vectors * sqrt(resolution)
  eigenvalues <- eig_Ce$values / resolution
  basis_tensors <- lapply(seq_len(resolution), function(k) {
    hkk <- as.vector(basis[, k] %o% basis[, k])
    hkk %o% hkk
  })

  alpha <- lapply(seq_len(p), function(i) {
    result <- matrix(0, nrow = resolution^2, ncol = resolution^2)
    for (k in seq_along(alpha_scale)) {
      result <- result + alpha_scale[k] * basis_tensors[[k]]
    }
    ap_amt[i] * result
  })
  if (p == 1) alpha <- alpha[[1]]
  delta <- C_eps

  errorsL2 <- errorsAbs <- errorsRel <- matrix(
    NA_real_, nrow = length(Ns), ncol = 2,
    dimnames = list(as.character(Ns), c("diag_new", "diag_new5"))
  )
  errorsL2_list <- errorsAbs_list <- errorsRel_list <- vector("list", length(Ns))

  build_alpha_operator <- function(alpha_values, fitted_basis, lag_index = 1) {
    fitted_alpha <- if (is.list(alpha_values)) {
      alpha_values[[lag_index]]
    } else if (is.null(dim(alpha_values))) {
      alpha_values
    } else {
      alpha_values[, lag_index]
    }
    result <- matrix(0, nrow = resolution^2, ncol = resolution^2)
    for (k in seq_along(fitted_alpha)) {
      result <- result + fitted_alpha[k] * basis_tensors[[k]]
    }
    result
  }

  alpha_error <- function(true_alpha, fitted_alpha, fitted_basis, fit_p) {
    true_list <- if (is.list(true_alpha)) true_alpha else list(true_alpha)
    fit_list <- lapply(seq_len(fit_p), function(i) {
      build_alpha_operator(fitted_alpha, fitted_basis, i)
    })
    n_compare <- max(length(true_list), length(fit_list))
    l2 <- abs_err <- denom <- vector("list", n_compare)

    for (i in seq_len(n_compare)) {
      truth <- if (i <= length(true_list)) true_list[[i]] else
        matrix(0, nrow = resolution^2, ncol = resolution^2)
      estimate <- if (i <= length(fit_list)) fit_list[[i]] else
        matrix(0, nrow = resolution^2, ncol = resolution^2)
      l2[[i]] <- pmin((truth - estimate)^2, (truth + estimate)^2)
      abs_err[[i]] <- pmin(abs(truth - estimate), abs(truth + estimate))
      denom[[i]] <- abs(truth)
    }

    list(
      L2 = fChange:::dot_integrate(Reduce("+", l2)),
      Abs = fChange:::dot_integrate(Reduce("+", abs_err)),
      Rel = fChange:::dot_integrate(Reduce("+", abs_err)) /
        fChange:::dot_integrate(Reduce("+", denom))
    )
  }

  for (idx in seq_along(Ns)) {
    N_current <- Ns[idx]
    if (!silent) cat("\nN=", N_current, ": ", sep = "")

    run_simulation <- function(sim) {
      simulated <- generateFARCH_use(
        N = N_current, resolution = resolution, alpha = alpha,
        delta = delta, basis = basis, burnin = burnin, p = p,
        eigs = eigenvalues[seq_along(alpha_scale)], C_eps = C_eps
      )

      diag_mod <- .fit_diag_model_new(
        X = simulated$data, C_eps = C_eps, basis = basis,
        eigenvalues = eigenvalues, w = w, p = 1
      )
      diag_mod5 <- .fit_diag_model_new(
        X = simulated$data, C_eps = C_eps, basis = basis,
        eigenvalues = eigenvalues, w = w, p = 5
      )

      if (parameter == "alpha") {
        err1 <- alpha_error(alpha, diag_mod$alpha, basis, fit_p = 1)
        err5 <- alpha_error(alpha, diag_mod5$alpha, basis, fit_p = 5)
      } else {
        err1_values <- pmin(
          (delta - diag_mod$delta)^2,
          (delta + diag_mod$delta)^2
        )
        err1_abs <- pmin(
          abs(delta - diag_mod$delta),
          abs(delta + diag_mod$delta)
        )
        err5_values <- pmin(
          (delta - diag_mod5$delta)^2,
          (delta + diag_mod5$delta)^2
        )
        err5_abs <- pmin(
          abs(delta - diag_mod5$delta),
          abs(delta + diag_mod5$delta)
        )
        delta_abs <- abs(delta)
        err1 <- list(
          L2 = fChange:::dot_integrate(err1_values),
          Abs = fChange:::dot_integrate(err1_abs),
          Rel = fChange:::dot_integrate(err1_abs) /
            fChange:::dot_integrate(delta_abs)
        )
        err5 <- list(
          L2 = fChange:::dot_integrate(err5_values),
          Abs = fChange:::dot_integrate(err5_abs),
          Rel = fChange:::dot_integrate(err5_abs) /
            fChange:::dot_integrate(delta_abs)
        )
      }

      c(L2_1 = err1$L2, L2_5 = err5$L2,
        Abs_1 = err1$Abs, Abs_5 = err5$Abs,
        Rel_1 = err1$Rel, Rel_5 = err5$Rel)
    }

    sim_results <- if (parallel) {
      future.apply::future_lapply(seq_len(nSims), run_simulation,
                                  future.seed = TRUE)
    } else {
      lapply(seq_len(nSims), run_simulation)
    }
    sim_results <- do.call(rbind, sim_results)

    errorsL2_list[[idx]] <- sim_results[, c("L2_1", "L2_5"), drop = FALSE]
    errorsAbs_list[[idx]] <- sim_results[, c("Abs_1", "Abs_5"), drop = FALSE]
    errorsRel_list[[idx]] <- sim_results[, c("Rel_1", "Rel_5"), drop = FALSE]
    colnames(errorsL2_list[[idx]]) <- colnames(errorsAbs_list[[idx]]) <-
      colnames(errorsRel_list[[idx]]) <- c("diag_new", "diag_new5")

    errorsL2[idx, ] <- colMeans(errorsL2_list[[idx]])
    errorsAbs[idx, ] <- colMeans(errorsAbs_list[[idx]])
    errorsRel[idx, ] <- colMeans(errorsRel_list[[idx]])

    if (!is.null(checkpoint_dir)) {
      saveRDS(
        list(parameter = parameter, p = p,
             completed_Ns = Ns[seq_len(idx)],
             L2 = errorsL2[seq_len(idx), , drop = FALSE],
             Abs = errorsAbs[seq_len(idx), , drop = FALSE],
             Rel = errorsRel[seq_len(idx), , drop = FALSE],
             lists = list(L2 = errorsL2_list[seq_len(idx)],
                          Abs = errorsAbs_list[seq_len(idx)],
                          Rel = errorsRel_list[seq_len(idx)])),
        file.path(checkpoint_dir,
                  paste0("convergence_", parameter, "_N", N_current, ".rds"))
      )
    }
  }

  list(
    parameter = parameter,
    p = p,
    L2 = errorsL2,
    Abs = errorsAbs,
    Rel = errorsRel,
    lists = list(L2 = errorsL2_list, Abs = errorsAbs_list,
                 Rel = errorsRel_list)
  )
}
