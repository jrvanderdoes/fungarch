#' Plot a functional residual surface
#'
#' Creates a three-dimensional Plotly line plot for a functional-data object.
#' The plotting function returns the widget; saving an image is left to the
#' calling script so that output paths and Python/Kaleido configuration remain
#' outside the package function.
#'
#' @param X A functional-data object with `data`, `fparam`, and `labels`
#'   components.
#' @param plot_title Optional plot title.
#' @param eye Plotly camera position.
#' @param aspectratio Plotly scene aspect ratio.
#' @param showticklabels Whether x-axis tick labels should be shown.
#' @return A Plotly `htmlwidget`.
#' @export
plot_residual_surface <- function(
    X, plot_title = NULL,
    eye = list(x = -0.5, y = -1.5, z = 0.25),
    aspectratio = list(x = 1.1, y = 0.5, z = 0.5),
    showticklabels = TRUE) {
  if (is.null(X$data) || is.null(X$fparam) || is.null(X$labels)) {
    stop("X must contain data, fparam, and labels components.")
  }
  if (ncol(X$data) != length(X$labels)) {
    stop("X$data and X$labels must have matching numbers of curves.")
  }

  plot_data <- do.call(rbind, lapply(seq_len(ncol(X$data)), function(i) {
    data.frame(
      resolution = X$fparam,
      FDRep = i,
      Value = X$data[, i]
    )
  }))

  plot_colors <- RColorBrewer::brewer.pal(11, "Spectral")
  plot_colors[6] <- "yellow"
  scene <- list(
    camera = list(eye = eye),
    aspectmode = "manual",
    aspectratio = aspectratio
  )

  plotly::plot_ly(
    plot_data,
    x = ~FDRep,
    y = ~resolution,
    z = ~Value,
    type = "scatter3d",
    mode = "lines",
    color = ~as.factor(FDRep),
    colors = plot_colors
  ) |>
    plotly::layout(scene = list(
      yaxis = list(
        title = "",
        showticklabels = FALSE,
        ticktext = round(seq(min(X$fparam), max(X$fparam), length.out = 5), 3),
        tickvals = fChange:::.select_n(
          vals = seq(min(X$fparam), max(X$fparam), length.out = 5), n = 5
        )
      ),
      xaxis = list(
        title = "",
        tickfont = list(size = 15),
        showticklabels = showticklabels,
        ticktext = fChange:::.select_n(
          vals = format(as.Date(X$labels), "%b %y"), n = 6
        ),
        tickvals = fChange:::.select_n(
          vals = seq_along(X$labels), n = 6
        )
      ),
      zaxis = list(title = "", showticklabels = FALSE)
    )) |>
    plotly::layout(title = plot_title, scene = scene) |>
    plotly::layout(showlegend = FALSE)
}
