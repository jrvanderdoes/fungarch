#' Forecast functional curves and prediction envelopes
#'
#' Fits the diagonal and pointwise functional ARCH models over a rolling
#' training period, computes lower and upper forecast bounds, and returns a
#' `ggplot2` comparison plot together with the calculated results.
#'
#' @param data Numeric matrix or data frame. Rows are discretized functional
#'   observations and columns are ordered curves.
#' @param N_fit Number of curves used before the first forecast.
#' @param forecast Number of forecast curves to evaluate.
#' @param w Smoothing/regularization weight passed to the model-fitting
#'   functions.
#' @param w1 Retained for compatibility; currently unused.
#' @param alp Tail probability used for the historical and model-based
#'   prediction bounds.
#' @return A list with components `estimates`, `plt`, `lower`, and `upper`.
#' @export
forecast_envelope <- function(data, N_fit = 450, forecast = 5,
                              w = 0.05, w1 = 0.05, alp = 0.05) {
  resolution <- nrow(data)
  C_eps <- .generate_ceps(resolution)
  eig_Ce <- eigen(C_eps)
  basis <- eig_Ce$vectors * sqrt(resolution)
  eigenvalues <- eig_Ce$values / resolution

  X2 <- array(apply(data, 2, tcrossprod),
              dim = c(nrow(data), nrow(data), ncol(data)))
  estimates <- estimates1 <- array(dim = c(resolution, 5, forecast))
  lower <- upper <- matrix(nrow = resolution, ncol = forecast)

  cat(paste0("(", forecast, ")"))
  # Kept as in the paper analysis: the first iteration uses the N_fit-th
  # curve as the last training curve and the next curve as the forecast.
  for (n in (N_fit - 1) + 1:forecast) {
    cat(paste0(n - N_fit + 1, ", "))
    X2_prev <- X2[, , n]
    X2_pred <- X2[, , n + 1]

    sds <- apply(data[, 1:n], MARGIN = 1, sd)
    means <- apply(data[, 1:n], MARGIN = 1, mean)

    diag_mod <- .fit_diag_model_new(
      X = data[, 1:n], C_eps = C_eps, basis = basis,
      eigenvalues = eigenvalues, w = w, p = 1
    )
    diag_pred <- .pred_diag_new(
      resolution, alphas = diag_mod$alpha, basis = basis,
      delta = diag_mod$delta, data = data[, 1:n], C_eps = C_eps
    )

    diag_mod5 <- .fit_diag_model_new(
      X = data[, 1:n], C_eps = C_eps, basis = basis,
      eigenvalues = eigenvalues, w = w, p = 5
    )
    diag_pred5 <- .pred_diag_new(
      resolution, alphas = diag_mod5$alpha, basis = basis,
      delta = diag_mod5$delta, data = data[, 1:n], C_eps = C_eps
    )

    pw_mod <- .fit_pw(X = data[, 1:n], w = w)
    pw_pred1 <- .pred_pw_simple(
      X = data[, 1:n], pw_mod$alpha, pw_mod$delta, C_eps
    )

    hist_pred <- apply(data[, 1:n], MARGIN = 1, quantile, probs = alp)
    mean_pred <- rowMeans(data[, 1:n])

    estimates[, , n - N_fit + 1] <- as.matrix(data.frame(
      data[, n + 1],
      sqrt(pmax(0, diag(diag_pred))) * qnorm(alp, sd = sds),
      sqrt(pmax(0, diag(diag_pred5))) * qnorm(alp, sd = sds),
      pw_pred1 * qnorm(alp, sd = sds),
      hist_pred
    ))
    estimates1[, , n - N_fit + 1] <- as.matrix(data.frame(
      data[, n + 1],
      sqrt(pmax(0, diag(diag_pred))) * qnorm(alp),
      sqrt(pmax(0, diag(diag_pred5))) * qnorm(alp),
      pw_pred1 * qnorm(alp),
      hist_pred
    ))

    lower[, n - N_fit + 1] <- sqrt(pmax(0, diag(diag_pred))) *
      qnorm(alp / 2, sd = sds)
    upper[, n - N_fit + 1] <- sqrt(pmax(0, diag(diag_pred))) *
      qnorm(1 - alp / 2, sd = sds)
  }

  tmp <- as.data.frame(estimates1[, 1, ])
  colnames(tmp) <- paste0("D", seq_len(forecast))
  tmp1 <- tmp |>
    dplyr::mutate(t = seq_len(resolution)) |>
    tidyr::pivot_longer(cols = paste0("D", seq_len(forecast)))
  tmp1 <- tmp1[order(tmp1$name, tmp1$t), , drop = FALSE]

  # The forecast labels are numbered after the fitting sample.
  x_labels <- paste0(N_fit + seq_len(forecast))
  if (forecast > 10) x_labels[seq(2, length(x_labels), by = 2)] <- ""

  cols <- scales::hue_pal()(4)
  plt <- ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(
      x = seq_len(nrow(tmp1)), y = c(estimates[, 1, ]),
      color = "True", group = rep(seq_len(forecast), each = resolution)
    ), linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(
      x = seq_len(nrow(tmp1)), y = c(estimates[, 2, ]),
      color = "CCC-op-ARCH(1)", group = rep(seq_len(forecast), each = resolution)
    )) +
    ggplot2::geom_line(ggplot2::aes(
      x = seq_len(nrow(tmp1)), y = c(estimates[, 3, ]),
      color = "CCC-op-ARCH(5)", group = rep(seq_len(forecast), each = resolution)
    )) +
    ggplot2::geom_line(ggplot2::aes(
      x = seq_len(nrow(tmp1)), y = c(estimates[, 4, ]),
      color = "pw-ARCH(1)", group = rep(seq_len(forecast), each = resolution)
    )) +
    ggplot2::geom_line(ggplot2::aes(
      x = seq_len(nrow(tmp1)), y = c(estimates[, 5, ]),
      color = "Historical", group = rep(seq_len(forecast), each = resolution)
    )) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = 0.5 + resolution * seq_len(forecast - 1)),
      linetype = "dashed"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = c(0.1, 0.15),
      axis.title = ggplot2::element_text(size = 20),
      axis.text = ggplot2::element_text(size = 16),
      legend.text = ggplot2::element_text(size = 16)
    ) +
    ggplot2::ylab("") +
    ggplot2::xlab("Day") +
    ggplot2::labs(color = NULL) +
    ggplot2::scale_x_continuous(
      breaks = resolution / 2 + resolution * 0:(forecast - 1),
      labels = x_labels
    ) +
    ggplot2::scale_color_manual(values = c(cols, "black"))

  list(estimates = estimates, plt = plt, lower = lower, upper = upper)
}

#' Plot observed curves and forecast envelopes
#'
#' @param X A functional-data object with `data` and `fparam` components.
#' @param lower Matrix of lower forecast bounds.
#' @param upper Matrix of upper forecast bounds.
#' @param changes Optional curve indices at which observed-data groups change.
#' @param plot_title Plot title.
#' @param val_axis_title Value-axis title.
#' @param res_axis_title Resolution-axis title.
#' @param FD_axis_title Functional-observation-axis title.
#' @param eye Plotly camera position.
#' @param aspectratio Optional Plotly scene aspect ratio.
#' @param showticklabels Whether axis tick labels should be shown.
#' @param fixme Whether to use the interpolated color palette for a single
#'   observed-data group.
#' @param ... Additional arguments reserved for compatibility.
#' @return A Plotly `htmlwidget`.
#' @export
plot_forecast_envelope <- function(
    X, lower, upper, changes = NULL, plot_title = X$name,
    val_axis_title = "Value", res_axis_title = "Resolution",
    FD_axis_title = "Observations",
    eye = list(x = 0.75, y = -1.5, z = 0.75), aspectratio = NULL,
    showticklabels = TRUE, fixme = TRUE, ...) {
  pred_n <- max(ncol(lower), ncol(upper))
  changes <- unique(c(changes, ncol(X) - pred_n))
  plotData <- plotData_forecast <- plotData_lower <- plotData_upper <- data.frame()

  if (length(changes) > 1) {
    for (j in 1:min(changes)) {
      plotData <- rbind(plotData, data.frame(
        resolution = X$fparam, FDRep = j, Color = 1, Value = X$data[, j]
      ))
    }
    for (i in 2:length(changes)) {
      for (j in (changes[i - 1] + 1):changes[i]) {
        plotData <- rbind(plotData, data.frame(
          resolution = X$fparam, FDRep = j, Color = i, Value = X$data[, j]
        ))
      }
    }
  } else {
    for (j in 1:min(changes)) {
      plotData <- rbind(plotData, data.frame(
        resolution = X$fparam, FDRep = j, Color = j, Value = X$data[, j]
      ))
    }
  }

  for (j in (max(changes) + 1):ncol(X)) {
    plotData_forecast <- rbind(plotData_forecast, data.frame(
      resolution = X$fparam, FDRep = j, Color = max(changes) + 1,
      Value = X$data[, j]
    ))
    plotData_lower <- rbind(plotData_lower, data.frame(
      resolution = X$fparam, FDRep = j, Color = max(changes) + 1,
      Value = lower[, j - max(changes)]
    ))
    plotData_upper <- rbind(plotData_upper, data.frame(
      resolution = X$fparam, FDRep = j, Color = max(changes) + 1,
      Value = upper[, j - max(changes)]
    ))
  }

  if (length(changes) > 1) {
    plot_colors <- RColorBrewer::brewer.pal(
      min(9, max(3, length(changes))), "Set1"
    )[1:min(9, length(changes))]
    if (length(changes) - 1 > 9) {
      plot_colors <- rep(
        plot_colors, ceiling((length(changes) + 1) / 9)
      )[1:length(changes)]
    }
    plot_colors <- c(plot_colors, "black")
    names(plot_colors) <- c(seq_len(length(changes)), max(changes) + 1)
  } else {
    plot_colors <- RColorBrewer::brewer.pal(11, "Spectral")
    plot_colors[6] <- "yellow"
    tmp <- if (fixme) {
      grDevices::colorRampPalette(plot_colors)(changes)
    } else {
      "black"
    }
    plot_colors <- c(tmp, "black", scales::hue_pal()(4)[2])
    names(plot_colors) <- seq_len(changes + 2)
  }

  if (is.null(aspectratio)) aspectratio <- list(x = 0.8, y = 0.5, z = 0.5)
  scene <- list(camera = list(eye = eye), aspectmode = "manual",
                aspectratio = aspectratio)

  plotly::plot_ly(
    rbind(plotData, plotData_forecast), x = ~as.factor(FDRep),
    y = ~resolution, z = ~Value, type = "scatter3d", mode = "lines",
    split = ~as.factor(FDRep), color = ~as.factor(Color), colors = plot_colors
  ) |>
    plotly::add_lines(
      data = plotData_lower, x = ~as.factor(FDRep), y = ~resolution,
      z = ~Value, type = "scatter3d", mode = "lines",
      split = ~as.factor(FDRep), color = as.character(length(plot_colors)),
      colors = plot_colors
    ) |>
    plotly::add_lines(
      data = plotData_upper, x = ~as.factor(FDRep), y = ~resolution,
      z = ~Value, type = "scatter3d", mode = "lines",
      split = ~as.factor(FDRep), color = as.character(length(plot_colors)),
      colors = plot_colors
    ) |>
    plotly::layout(
      scene = list(
        yaxis = list(title = res_axis_title, showticklabels = showticklabels),
        xaxis = list(title = FD_axis_title, showticklabels = showticklabels),
        zaxis = list(title = val_axis_title, showticklabels = showticklabels)
      )
    ) |>
    plotly::layout(title = plot_title, scene = scene) |>
    plotly::layout(showlegend = FALSE)
}
