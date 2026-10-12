# Tables comparing one-step pointwise forecast violations for SPY data.
#
# This combines the former table_spy_1820.R and table_spy_2224.R scripts.
# The original scripts are no longer part of the package; this is the maintained
# paper workflow.

# if (!requireNamespace("fungarch", quietly = TRUE)) {
#   stop(
#     "Install fungarch before running this script with devtools::install('.').",
#     call. = FALSE
#   )
# }
dir.create("./paper/figs", recursive = TRUE, showWarnings = FALSE)

prepare_spy_table_data <- function(x) {
  ocidr_data <- fChange::ocidr(x)$data
  colnames(ocidr_data) <- colnames(x)
  ocidr_data
}

# Create one independent L'Ecuyer-CMRG stream per forecast position.  Using
# the forecast position, rather than the position within a chunk, makes
# results reproducible when the work is split across computers.
make_spy_rng_seeds <- function(indices, seed = 12345L) {
  indices <- as.integer(indices)
  if (length(indices) == 0L) return(list())

  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(as.integer(seed))
  all_seeds <- vector("list", max(indices))
  all_seeds[[1L]] <- .Random.seed
  if (length(all_seeds) > 1L) {
    for (j in 2:length(all_seeds)) {
      all_seeds[[j]] <- parallel::nextRNGStream(all_seeds[[j - 1L]])
    }
  }
  all_seeds[indices]
}

run_spy_violation_table <- function(
    data,
    training_end,
    w = 0.05,
    alp = 0.05,
    silent = FALSE) {
  if (!is.matrix(data) || !is.numeric(data)) {
    stop("data must be a numeric matrix.", call. = FALSE)
  }
  if (length(training_end) != 1L || !inherits(training_end, "Date")) {
    stop("training_end must be one Date.", call. = FALSE)
  }
  if (length(alp) != 1L || !is.numeric(alp) ||
      !is.finite(alp) || alp <= 0 || alp >= 1) {
    stop("alp must be one number strictly between 0 and 1.", call. = FALSE)
  }

  dates <- as.Date(colnames(data), format = "%m-%d-%Y")
  if (anyNA(dates)) {
    stop("All data columns must have names in MM-DD-YYYY format.",
         call. = FALSE)
  }
  N_fit <- sum(dates < training_end)
  if (N_fit < 3L || ncol(data) - N_fit - 1L < 1L) {
    stop("training_end must leave at least one test day.", call. = FALSE)
  }

  resolution <- nrow(data)
  C_eps <- fungarch:::.generate_ceps(resolution)
  eig_Ce <- eigen(C_eps)
  basis <- eig_Ce$vectors * sqrt(resolution)
  eigenvalues <- eig_Ce$values / resolution

  n_forecasts <- ncol(data) - N_fit - 1L
  violations <- matrix(
    NA_real_, nrow = n_forecasts, ncol = 4,
    dimnames = list(NULL, c("diag", "diag5", "pointwise", "historical"))
  )

  for (i in seq_len(n_forecasts)) {
    # Match the original scripts' first forecast origin.
    n <- N_fit + i
    if (!silent) cat(i, ",", sep = "")

    X_train <- data[, seq_len(n), drop = FALSE]
    X_next <- data[, n + 1L]
    diag_mod <- fungarch:::.fit_diag_model_new(
      X = X_train, C_eps = C_eps, basis = basis,
      eigenvalues = eigenvalues, w = w, p = 1
    )
    diag_pred <- fungarch:::.pred_diag_new(
      resolution, alphas = diag_mod$alpha, basis = basis,
      delta = diag_mod$delta, data = X_train, C_eps = C_eps
    )

    diag_mod5 <- fungarch:::.fit_diag_model_new(
      X = X_train, C_eps = C_eps, basis = basis,
      eigenvalues = eigenvalues, w = w, p = 5
    )
    diag_pred5 <- fungarch:::.pred_diag_new(
      resolution, alphas = diag_mod5$alpha, basis = basis,
      delta = diag_mod5$delta, data = X_train, C_eps = C_eps
    )

    pw_mod <- fungarch:::.fit_pw(X = X_train, w = w)
    pw_pred <- fungarch:::.pred_pw_simple(
      X = X_train, pw_mod$alpha, pw_mod$delta, C_eps
    )

    historical_pred <- apply(X_train, 1, quantile, probs = alp)
    q_threshold <- qnorm(alp)

    predictions <- cbind(
      sqrt(pmax(0, diag(diag_pred))) * q_threshold,
      sqrt(pmax(0, diag(diag_pred5))) * q_threshold,
      pw_pred * q_threshold,
      historical_pred
    )

    violations[i, ] <- colMeans(
      X_next <= predictions,
      na.rm = TRUE
    )
  }

  if (!silent) cat("\n")

  list(
    violations = violations,
    summary = colMeans(violations, na.rm = TRUE),
    N_fit = N_fit,
    alpha = alp,
    w = w,
    training_end = training_end
  )
}

run_spy_table <- function(
    data,
    training_end,
    w,
    alpha_levels = c(0.05, 0.01),
    silent = FALSE,
    forecast_indices = NULL,
    parallel = FALSE,
    ncores = min(4L, max(1L, parallel::detectCores(logical = TRUE) - 1L)),
    checkpoint_file = NULL,
    seed = 12345L) {
  if (length(alpha_levels) < 1L || any(!is.finite(alpha_levels)) ||
      any(alpha_levels <= 0 | alpha_levels >= 1)) {
    stop("alpha_levels must contain values strictly between 0 and 1.",
         call. = FALSE)
  }

  dates <- as.Date(colnames(data), format = "%m-%d-%Y")
  if (anyNA(dates)) {
    stop("All data columns must have names in MM-DD-YYYY format.",
         call. = FALSE)
  }
  N_fit <- sum(dates < training_end)
  n_forecasts <- ncol(data) - N_fit - 1L
  if (N_fit < 3L || n_forecasts < 1L) {
    stop("training_end must leave at least one test day.", call. = FALSE)
  }
  if (is.null(forecast_indices)) {
    forecast_indices <- seq_len(n_forecasts)
  }
  if (!is.numeric(forecast_indices) ||
      any(forecast_indices != as.integer(forecast_indices)) ||
      any(forecast_indices < 1 | forecast_indices > n_forecasts)) {
    stop("forecast_indices must be valid forecast positions.", call. = FALSE)
  }
  forecast_indices <- as.integer(forecast_indices)
  if (!is.logical(parallel) || length(parallel) != 1L) {
    stop("parallel must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(ncores) || length(ncores) != 1L ||
      ncores != as.integer(ncores) || ncores < 1) {
    stop("ncores must be one positive integer.", call. = FALSE)
  }
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("seed must be one finite integer.", call. = FALSE)
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

  rng_seeds <- make_spy_rng_seeds(forecast_indices, seed = seed)

  resolution <- nrow(data)
  C_eps <- fungarch:::.generate_ceps(resolution)
  eig_Ce <- eigen(C_eps)
  basis <- eig_Ce$vectors * sqrt(resolution)
  eigenvalues <- eig_Ce$values / resolution

  violation_names <- c("diag", "diag5", "pointwise", "historical")
  fit_one_forecast <- function(i) {
    # Match the original scripts' first forecast origin.
    n <- N_fit + i

    X_train <- data[, seq_len(n), drop = FALSE]
    X_next <- data[, n + 1L]
    # These expensive fits and predictions are shared by both alpha levels.
    diag_mod <- fungarch:::.fit_diag_model_new(
      X = X_train, C_eps = C_eps, basis = basis,
      eigenvalues = eigenvalues, w = w, p = 1
    )
    diag_pred <- fungarch:::.pred_diag_new(
      resolution, alphas = diag_mod$alpha, basis = basis,
      delta = diag_mod$delta, data = X_train, C_eps = C_eps
    )

    diag_mod5 <- fungarch:::.fit_diag_model_new(
      X = X_train, C_eps = C_eps, basis = basis,
      eigenvalues = eigenvalues, w = w, p = 5
    )
    diag_pred5 <- fungarch:::.pred_diag_new(
      resolution, alphas = diag_mod5$alpha, basis = basis,
      delta = diag_mod5$delta, data = X_train, C_eps = C_eps
    )

    pw_mod <- fungarch:::.fit_pw(X = X_train, w = w)
    pw_pred <- fungarch:::.pred_pw_simple(
      X = X_train, pw_mod$alpha, pw_mod$delta, C_eps
    )

    rows <- lapply(seq_along(alpha_levels), function(a) {
      alp <- alpha_levels[a]
      q_threshold <- qnorm(alp)
      historical_pred <- apply(X_train, 1, quantile, probs = alp)

      predictions <- cbind(
        sqrt(pmax(0, diag(diag_pred))) * q_threshold,
        sqrt(pmax(0, diag(diag_pred5))) * q_threshold,
        pw_pred * q_threshold,
        historical_pred
      )

      colMeans(
        X_next <= predictions,
        na.rm = TRUE
      )
    })
    names(rows) <- paste0("alpha_", alpha_levels)
    rows
  }

  if (parallel) {
    forecast_results <- future.apply::future_lapply(
      forecast_indices, fit_one_forecast, future.seed = rng_seeds
    )
  } else {
    old_had_seed <- exists(".Random.seed", envir = .GlobalEnv,
                           inherits = FALSE)
    old_seed <- if (old_had_seed) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (old_had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv,
                        inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    forecast_results <- lapply(seq_along(forecast_indices), function(j) {
      assign(".Random.seed", rng_seeds[[j]], envir = .GlobalEnv)
      fit_one_forecast(forecast_indices[[j]])
    })
  }

  violations <- lapply(seq_along(alpha_levels), function(a) {
    result <- do.call(rbind, lapply(forecast_results, `[[`, a))
    colnames(result) <- violation_names
    rownames(result) <- as.character(forecast_indices)
    result
  })
  names(violations) <- paste0("alpha_", alpha_levels)

  if (!silent) cat("\n")

  results <- lapply(seq_along(alpha_levels), function(a) {
    list(
      violations = violations[[a]],
      summary = colMeans(violations[[a]], na.rm = TRUE),
      forecast_indices = forecast_indices,
      n_forecasts = n_forecasts,
      N_fit = N_fit,
      alpha = alpha_levels[a],
      w = w,
      training_end = training_end,
      seed = seed
    )
  })
  names(results) <- paste0("alpha_", alpha_levels)
  if (!is.null(checkpoint_file)) {
    dir.create(dirname(checkpoint_file), recursive = TRUE,
               showWarnings = FALSE)
    saveRDS(results, checkpoint_file)
  }
  results
}

combine_spy_table_results <- function(...) {
  chunks <- list(...)
  if (length(chunks) < 1L) {
    stop("Supply at least one result chunk.", call. = FALSE)
  }
  alpha_names <- names(chunks[[1]])
  if (is.null(alpha_names) ||
      any(!vapply(chunks, function(x) identical(names(x), alpha_names),
                  logical(1)))) {
    stop("All chunks must contain the same alpha levels.", call. = FALSE)
  }

  combined <- lapply(alpha_names, function(alpha_name) {
    pieces <- lapply(chunks, function(chunk) {
      chunk[[alpha_name]]$violations
    })
    violations <- do.call(rbind, pieces)
    violations <- violations[
      order(as.integer(rownames(violations))), , drop = FALSE
    ]
    result <- chunks[[1]][[alpha_name]]
    result$violations <- violations
    result$forecast_indices <- as.integer(rownames(violations))
    result$summary <- colMeans(violations, na.rm = TRUE)
    result
  })
  names(combined) <- alpha_names
  combined
}

data("SPY_18_20", package = "fungarch")
data("SPY_22_24", package = "fungarch")

spy_table_results <- list(
  SPY_18_20 = run_spy_table(
    prepare_spy_table_data(SPY_18_20),
    training_end = as.Date("2020-01-01"),
    w = 0.05,
    parallel = TRUE,
    ncores = min(4L, max(1L, parallel::detectCores(logical = TRUE) - 1L)),
    checkpoint_file = "./paper/figs/table_spy_18_20.rds"
  ),
  SPY_22_24 = run_spy_table(
    prepare_spy_table_data(SPY_22_24),
    training_end = as.Date("2024-01-01"),
    w = 0.10,
    parallel = TRUE,
    ncores = min(4L, max(1L, parallel::detectCores(logical = TRUE) - 1L)),
    checkpoint_file = "./paper/figs/table_spy_22_24.rds"
  )
)

print(lapply(spy_table_results, function(x) {
  lapply(x, `[[`, "summary")
}))

saveRDS(spy_table_results, "./paper/figs/table_spy_results.rds")

# To split SPY_18_20 across two computers, use the same data and package
# versions on both machines, and divide the forecast positions. Because the
# calculation is deterministic, no random seed is needed.
#
# spy_18_20_data <- prepare_spy_table_data(SPY_18_20)
# n_fit_18_20 <- sum(
#   as.Date(colnames(spy_18_20_data), "%m-%d-%Y") < as.Date("2020-01-01")
# )
# n_forecasts_18_20 <- ncol(spy_18_20_data) - n_fit_18_20
# split_18_20 <- floor(n_forecasts_18_20 / 2)
#
# part1 <- run_spy_table(
#   spy_18_20_data, as.Date("2020-01-01"), w = 0.05,
#   forecast_indices = seq_len(split_18_20),
#   parallel = TRUE, ncores = 6,
#   checkpoint_file = "./paper/figs/table_spy_18_20_part1.rds"
# )
# part2 <- run_spy_table(
#   spy_18_20_data, as.Date("2020-01-01"), w = 0.05,
#   forecast_indices = seq.int(split_18_20 + 1, n_forecasts_18_20),
#   parallel = TRUE, ncores = 6,
#   checkpoint_file = "./paper/figs/table_spy_18_20_part2.rds"
# )
# combined_18_20 <- combine_spy_table_results(part1, part2)
# print(lapply(combined_18_20, `[[`, "summary"))
#
# The same pattern works for SPY_22_24; change the dataset, training date,
# and `w` value:
#
# spy_22_24_data <- prepare_spy_table_data(SPY_22_24)
# n_fit_22_24 <- sum(
#   as.Date(colnames(spy_22_24_data), "%m-%d-%Y") < as.Date("2024-01-01")
# )
# n_forecasts_22_24 <- ncol(spy_22_24_data) - n_fit_22_24
# split_22_24 <- floor(n_forecasts_22_24 / 2)
#
# part1 <- run_spy_table(
#   spy_22_24_data, as.Date("2024-01-01"), w = 0.10,
#   forecast_indices = seq_len(split_22_24),
#   parallel = TRUE, ncores = 6,
#   checkpoint_file = "./paper/figs/table_spy_22_24_part1.rds"
# )
# part2 <- run_spy_table(
#   spy_22_24_data, as.Date("2024-01-01"), w = 0.10,
#   forecast_indices = seq.int(split_22_24 + 1, n_forecasts_22_24),
#   parallel = TRUE, ncores = 6,
#   checkpoint_file = "./paper/figs/table_spy_22_24_part2.rds"
# )
# combined_22_24 <- combine_spy_table_results(part1, part2)
# print(lapply(combined_22_24, `[[`, "summary"))
