devtools::load_all()
data("SPY_18_20", package = "fungarch")
dir.create("./paper/figs", recursive = TRUE, showWarnings = FALSE)

python_path <- Sys.getenv("RETICULATE_PYTHON")
if (nzchar(python_path)) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("RETICULATE_PYTHON is set, but the reticulate package is not installed.")
  }
  reticulate::use_python(python_path, required = TRUE)
}

data <- fChange::ocidr(SPY_18_20)$data
colnames(data) <- colnames(SPY_18_20)

res <- compute_white_noise_residuals(data = data)

for (lag in c(3, 10)) {
  print(fChange::portmanteau_tests(res$x_hat^2, test = "multi", lag = lag)$pvalue)
  print(fChange::portmanteau_tests(res$e_hat^2, test = "multi", lag = lag)$pvalue)
  print(fChange::portmanteau_tests(res$e_hat5^2, test = "multi", lag = lag)$pvalue)
}

X <- fChange::dfts(
  res$e_hat,
  labels = as.Date(colnames(data), "%m-%d-%Y")[-1],
  inc.warnings = FALSE
)
plotly::save_image(plot_residual_surface(X),
                   file = "./paper/figs/sp500_resids_1_full.png")

X <- fChange::dfts(
  res$e_hat5,
  labels = as.Date(colnames(data), "%m-%d-%Y")[-c(1:5)],
  inc.warnings = FALSE
)
plotly::save_image(plot_residual_surface(X),
                   file = "./paper/figs/sp500_resids_5.png")

png("./paper/figs/x_orig_wn.png", width = 800, height = 800)
print(fChange::sacf(data, lag.max = 10)$plot +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 40),
                 axis.title = ggplot2::element_text(size = 50)))
dev.off()

png("./paper/figs/x_hat_wn.png", width = 800, height = 800)
print(fChange::sacf(res$x_hat^2, lag.max = 10)$plot +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 40),
                 axis.title = ggplot2::element_text(size = 50)))
dev.off()

png("./paper/figs/e_hat_wn.png", width = 800, height = 800)
print(fChange::sacf(res$e_hat^2, lag.max = 10)$plot +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 40),
                 axis.title = ggplot2::element_text(size = 50)))
dev.off()

png("./paper/figs/e_hat5_wn.png", width = 800, height = 800)
print(fChange::sacf(res$e_hat5^2, lag.max = 10)$plot +
  ggplot2::theme(axis.text = ggplot2::element_text(size = 40),
                 axis.title = ggplot2::element_text(size = 50)))
dev.off()
