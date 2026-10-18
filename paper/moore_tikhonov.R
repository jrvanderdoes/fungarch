# Compare Moore--Penrose and Tikhonov diagonal estimators.
#
# This is the faster maintained Moore--Tikhonov workflow. It retains the rolling
# estimation design used for the figure, but omits covariance predictions and
# coverage calculations that are not used in the final relative-error plots.

devtools::load_all()

build_alpha_operator_moore_tikhonov <- function(
    alpha_values, basis = NULL, basis_tensors = NULL) {
  alpha_values <- as.numeric(alpha_values)
  if (is.null(basis) && is.null(basis_tensors)) {
    stop("Supply basis or basis_tensors.", call. = FALSE)
  }

  if (!is.null(basis)) {
    if (length(alpha_values) > ncol(basis)) {
      stop("alpha_values cannot exceed the number of basis columns.",
           call. = FALSE)
    }
    resolution <- nrow(basis)
    result <- matrix(0, nrow = resolution^2, ncol = resolution^2)
    for (k in seq_along(alpha_values)) {
      hkk <- as.vector(basis[, k] %o% basis[, k])
      result <- result + alpha_values[[k]] * (hkk %o% hkk)
    }
    return(result)
  }

  Reduce(
    `+`,
    Map(`*`, alpha_values[seq_along(basis_tensors)], basis_tensors),
    init = matrix(0, nrow = nrow(basis_tensors[[1]]),
                  ncol = ncol(basis_tensors[[1]]))
  )
}

relative_alpha_error_moore_tikhonov <- function(truth, estimate) {
  abs_error <- pmin(abs(truth - estimate), abs(truth + estimate))
  denominator <- abs(truth)
  fChange:::dot_integrate(abs_error) /
    fChange:::dot_integrate(denominator)
}

prepare_moore_tikhonov <- function(resolution, alpha_scale) {
  C_eps <- fungarch:::.generate_ceps(resolution)
  eig_Ce <- eigen(C_eps)
  basis <- eig_Ce$vectors * sqrt(resolution)
  eigenvalues <- eig_Ce$values / resolution
  alpha <- build_alpha_operator_moore_tikhonov(
    alpha_scale, basis = basis
  )

  list(
    C_eps = C_eps,
    basis = basis,
    eigenvalues = eigenvalues,
    alpha = alpha
  )
}

compare_models_moore_tikhonov <- function(
    N = 750,
    resolution = 50,
    alpha_scale = c(0.7, 0.7),
    N_fit = 500,
    w = 0.1,
    silent = FALSE,
    burnin = 100,
    setup = NULL) {
  if (N_fit < 2 || N_fit >= N) {
    stop("N_fit must be at least 2 and smaller than N.", call. = FALSE)
  }
  if (length(alpha_scale) > resolution) {
    stop("alpha_scale cannot contain more values than resolution.", call. = FALSE)
  }

  if (is.null(setup)) {
    setup <- prepare_moore_tikhonov(resolution, alpha_scale)
  }
  C_eps <- setup$C_eps
  basis <- setup$basis
  eigenvalues <- setup$eigenvalues
  alpha <- setup$alpha

  simulated <- fungarch:::generateFARCH_use(
    N = N,
    resolution = resolution,
    alpha = alpha,
    delta = C_eps,
    basis = basis,
    burnin = burnin,
    p = 1,
    eigs = eigenvalues[seq_along(alpha_scale)],
    C_eps = C_eps
  )

  n_forecasts <- N - N_fit
  relative_error <- matrix(
    NA_real_, nrow = n_forecasts, ncol = 2,
    dimnames = list(NULL, c("Moore-Penrose", "Tikhonov"))
  )

  for (i in seq_len(n_forecasts)) {
    forecast_index <- N_fit + i
    if (!silent) cat(i, ",", sep = "")

    X_fit <- simulated$data[, seq_len(forecast_index - 1), drop = FALSE]

    mp_model <- fungarch:::.fit_diag_model_new(
      X = X_fit,
      C_eps = C_eps,
      basis = basis,
      eigenvalues = eigenvalues,
      w = w,
      p = 1
    )
    tik_model <- fungarch:::.fit_diag_model_tik(
      X = X_fit,
      C_eps = C_eps,
      basis = basis,
      eigenvalues = eigenvalues,
      w = w,
      p = 1
    )

    mp_alpha <- build_alpha_operator_moore_tikhonov(
      mp_model$alpha, basis = basis
    )
    tik_alpha <- build_alpha_operator_moore_tikhonov(
      tik_model$alpha, basis = basis
    )

    relative_error[i, ] <- c(
      relative_alpha_error_moore_tikhonov(alpha, mp_alpha),
      relative_alpha_error_moore_tikhonov(alpha, tik_alpha)
    )
  }

  if (!silent) cat("\n")
  relative_error
}

run_one_moore_tikhonov_simulation <- function(
    simulation,
    N,
    resolution,
    alpha_scale,
    N_fit,
    w,
    setup = NULL,
    silent = TRUE) {
  if (is.null(setup)) {
    setup <- prepare_moore_tikhonov(resolution, alpha_scale)
  }
  errors <- compare_models_moore_tikhonov(
    N = N,
    resolution = resolution,
    alpha_scale = alpha_scale,
    N_fit = N_fit,
    w = w,
    silent = silent,
    setup = setup
  )
  colMeans(errors)
}

run_moore_tikhonov_study <- function(
    sizes = c(50, 100, 250, 500, 750, 1000),
    nSims = 50,
    resolution = 50,
    alpha_scale = c(0.7, 0.7),
    N_fit_offset = 10,
    w = 0.1,
    seed = 12345,
    silent = FALSE,
    parallel = FALSE,
    ncores = min(4L, max(1L, parallel::detectCores(logical = TRUE) - 1L)),
    checkpoint_file = NULL,
    resume = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  if (!is.logical(parallel) || length(parallel) != 1) {
    stop("parallel must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(ncores) || length(ncores) != 1 ||
      ncores != as.integer(ncores) || ncores < 1) {
    stop("ncores must be one positive integer.", call. = FALSE)
  }
  if (!is.logical(resume) || length(resume) != 1 || is.na(resume)) {
    stop("resume must be TRUE or FALSE.", call. = FALSE)
  }
  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop("parallel = TRUE requires future and future.apply.", call. = FALSE)
    }
    old_plan <- future::plan()
    future::plan(future::multisession, workers = as.integer(ncores))
    on.exit(future::plan(old_plan), add = TRUE)
  }

  result <- matrix(
    NA_real_, nrow = length(sizes), ncol = 2,
    dimnames = list(as.character(sizes), c("Moore-Penrose", "Tikhonov"))
  )

  # This setup depends only on the study settings, not on N or the
  # simulation replicate. Build it once and reuse it everywhere.
  setup <- prepare_moore_tikhonov(resolution, alpha_scale)

  # Checkpoints are written after each completed data length.  Restore only
  # compatible completed rows; the per-size seed below ensures that resumed
  # and uninterrupted runs generate the same simulations.
  if (resume && !is.null(checkpoint_file) &&
      file.exists(checkpoint_file)) {
    checkpoint <- readRDS(checkpoint_file)
    required_fields <- c("result", "completed_sizes", "nSims",
                         "resolution", "alpha_scale", "w")
    if (!all(required_fields %in% names(checkpoint))) {
      stop("Checkpoint is missing required study settings.", call. = FALSE)
    }
    same <- function(x, y) isTRUE(all.equal(x, y, check.attributes = FALSE))
    if (!same(checkpoint$nSims, nSims) ||
        !same(checkpoint$resolution, resolution) ||
        !same(checkpoint$alpha_scale, alpha_scale) ||
        !same(checkpoint$w, w)) {
      stop("Checkpoint settings do not match the requested study.",
           call. = FALSE)
    }

    completed_sizes <- as.numeric(checkpoint$completed_sizes)
    if (anyNA(completed_sizes) || any(!completed_sizes %in% sizes)) {
      stop("Checkpoint contains data lengths not present in `sizes`.",
           call. = FALSE)
    }
    if (nrow(checkpoint$result) != length(completed_sizes) ||
        ncol(checkpoint$result) != ncol(result)) {
      stop("Checkpoint result dimensions do not match completed sizes.",
           call. = FALSE)
    }

    result[match(completed_sizes, sizes), ] <- checkpoint$result
    if (!silent) {
      cat("Resuming from checkpoint; completed sizes: ",
          paste(completed_sizes, collapse = ", "), "\n", sep = "")
    }
  }

  for (i in seq_along(sizes)) {
    if (all(!is.na(result[i, ]))) {
      if (!silent) cat("Skipping completed ", sizes[i], "\n", sep = "")
      next
    }

    # Derive a stable stream from the sample size. This makes a given N
    # reproducible whether it is run alone or as part of the full study.
    seed_for_size <- if (is.null(seed)) {
      NULL
    } else {
      as.integer((abs(as.numeric(seed)) + as.numeric(sizes[i])) %%
                   2147483646)
    }
    if (!silent) cat("\nData length = ", sizes[i], "\n", sep = "")
    if (parallel) {
      simulation_errors <- future.apply::future_lapply(
        seq_len(nSims),
        run_one_moore_tikhonov_simulation,
        N = sizes[i],
        resolution = resolution,
        alpha_scale = alpha_scale,
        N_fit = sizes[i] - N_fit_offset,
        w = w,
        setup = setup,
        silent = TRUE,
        future.seed = if (is.null(seed_for_size)) FALSE else seed_for_size
      )
    } else {
      if (!is.null(seed_for_size)) set.seed(seed_for_size)
      simulation_errors <- lapply(
        seq_len(nSims),
        run_one_moore_tikhonov_simulation,
        N = sizes[i],
        resolution = resolution,
        alpha_scale = alpha_scale,
        N_fit = sizes[i] - N_fit_offset,
        w = w,
        setup = setup,
        silent = TRUE
      )
    }
    simulation_errors <- do.call(rbind, simulation_errors)

    result[i, ] <- colMeans(simulation_errors)

    if (!is.null(checkpoint_file)) {
      dir.create(dirname(checkpoint_file), recursive = TRUE,
                 showWarnings = FALSE)
      saveRDS(
        list(result = result[seq_len(i), , drop = FALSE],
             completed_sizes = sizes[seq_len(i)],
             nSims = nSims, resolution = resolution,
             alpha_scale = alpha_scale, w = w, seed = seed,
             sizes = sizes),
        checkpoint_file
      )
    }

    if (!silent) cat("Completed ", sizes[i], "\n", sep = "")
  }

  result
}

plot_moore_tikhonov <- function(low_dimensional, high_dimensional,
                        file = "./paper/figs/MPvT.png") {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  plot_data <- rbind(
    data.frame(
      DataLength = as.numeric(rownames(low_dimensional)),
      RelativeError = low_dimensional[, "Moore-Penrose"],
      Estimator = "Moore-Penrose",
      Setting = "Low-dimensional"
    ),
    data.frame(
      DataLength = as.numeric(rownames(low_dimensional)),
      RelativeError = low_dimensional[, "Tikhonov"],
      Estimator = "Tikhonov",
      Setting = "Low-dimensional"
    ),
    data.frame(
      DataLength = as.numeric(rownames(high_dimensional)),
      RelativeError = high_dimensional[, "Moore-Penrose"],
      Estimator = "Moore-Penrose",
      Setting = "High-dimensional"
    ),
    data.frame(
      DataLength = as.numeric(rownames(high_dimensional)),
      RelativeError = high_dimensional[, "Tikhonov"],
      Estimator = "Tikhonov",
      Setting = "High-dimensional"
    )
  )

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = DataLength, y = RelativeError,
                 color = Estimator, group = Estimator)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~Setting, nrow = 1) +
    ggplot2::scale_color_manual(
      values = c("Moore-Penrose" = "#F8766D", "Tikhonov" = "#2CB8C4")
    ) +
    ggplot2::labs(
      x = "Data Length",
      y = "Relative Absolute Error",
      color = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = c(0.82, 0.88),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank()
    )

  ggplot2::ggsave(file, plot = plot, width = 1600, height = 800,
                  units = "px", dpi = 100)
  invisible(plot)
}

# The two runs correspond to the left and right panels in the paper figure.
ncores <- min(4L, max(1L, parallel::detectCores(logical = TRUE) - 1L))
low_dimensional <- run_moore_tikhonov_study(
  alpha_scale = c(0.7, 0.7),
  parallel = TRUE,
  ncores = ncores,
  checkpoint_file = "./paper/figs/moore_tikhonov_low.rds"
)
high_dimensional <- run_moore_tikhonov_study(
  alpha_scale = 1 / (1:20)^2,
  parallel = TRUE,
  ncores = ncores,
  checkpoint_file = "./paper/figs/moore_tikhonov_high.rds"
)
plot_moore_tikhonov(low_dimensional, high_dimensional)

# Optional split across computers. Run one block on each machine, then combine
# the returned matrices with rbind() and sort by row names. The stable
# per-sample-size seed above makes the results match a single full run.
#
# low_part_1 <- run_moore_tikhonov_study(
#   sizes = c(50, 100, 250),
#   alpha_scale = c(0.7, 0.7),
#   parallel = TRUE,
#   ncores = 6,
#   checkpoint_file = "./paper/figs/moore_tikhonov_low_part1.rds"
# )
# low_part_2 <- run_moore_tikhonov_study(
#   sizes = c(500, 750, 1000),
#   alpha_scale = c(0.7, 0.7),
#   parallel = TRUE,
#   ncores = 6,
#   checkpoint_file = "./paper/figs/moore_tikhonov_low_part2.rds"
# )
# low_dimensional <- rbind(low_part_1, low_part_2)
# low_dimensional <- low_dimensional[
#   order(as.numeric(rownames(low_dimensional))), , drop = FALSE
# ]
