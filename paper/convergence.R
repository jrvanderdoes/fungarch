devtools::load_all()
dir.create("./paper/figs", recursive = TRUE, showWarnings = FALSE)

# One seed controls all scenarios reproducibly. The original scripts used
# separate seeds, so the individual curves will not be identical to every
# historical run unless the same seed and settings are restored there.
set.seed(123)

ncores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
checkpoint_dir <- "./paper/figs/convergence_checkpoints"

plot_convergence <- function(result, file) {
  plot_data <- as.data.frame(result$Rel)
  plot_data$n <- as.numeric(rownames(plot_data))
  plot_data <- tidyr::pivot_longer(
    plot_data,
    cols = c("diag_new", "diag_new5"),
    names_to = "model",
    values_to = "value"
  )

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = n, y = value, group = model, color = model)
  ) +
    ggplot2::geom_line(linewidth = 1.5) +
    ggplot2::theme_bw() +
    ggplot2::ylab("Relative Absolute Error") +
    ggplot2::xlab("Data Length") +
    ggplot2::scale_color_hue(
      labels = c("CCC-op-ARCH(1)", "CCC-op-ARCH(5)")
    ) +
    ggplot2::labs(color = NULL) +
    ggplot2::theme(
      legend.position = c(0.775, 0.85),
      axis.title = ggplot2::element_text(size = 35),
      axis.text = ggplot2::element_text(size = 30),
      legend.text = ggplot2::element_text(size = 30),
      plot.margin = ggplot2::margin(t = 10, b = 10, l = 10, r = 20)
    )

  ggplot2::ggsave(file, plot = plot, width = 800, height = 800,
                  units = "px", dpi = 100)
  invisible(plot)
}

results_A1 <- compare_convergence(
  parameter = "alpha",
  alpha_scale = c(0.7, 0.5),
  ap_amt = 1,
  p = 1,
  parallel = TRUE,
  ncores = ncores,
  checkpoint_dir = checkpoint_dir
)
print(results_A1$Rel)
plot_convergence(results_A1, "./paper/figs/convergence_A1.png")

results_A1_delta <- compare_convergence(
  parameter = "delta",
  alpha_scale = 0.6,
  ap_amt = 1,
  p = 1,
  parallel = TRUE,
  ncores = ncores,
  checkpoint_dir = checkpoint_dir
)
print(results_A1_delta$Rel)
plot_convergence(results_A1_delta, "./paper/figs/convergence_A1_delta.png")

results_A5 <- compare_convergence(
  parameter = "alpha",
  resolution = 20,
  alpha_scale = c(0.7, 0.5),
  ap_amt = c(0.1, 0.1, 0.1, 0.1, 0.3),
  p = 5,
  parallel = TRUE,
  ncores = ncores,
  checkpoint_dir = checkpoint_dir
)
print(results_A5$Rel)
plot_convergence(results_A5, "./paper/figs/convergence_A5.png")

saveRDS(
  list(A1 = results_A1, A1_delta = results_A1_delta, A5 = results_A5),
  "./paper/figs/convergence_results.rds"
)
