devtools::load_all()
data("SPY_22_24", package = "fungarch")
dir.create("./paper/figs", recursive = TRUE, showWarnings = FALSE)

data_22_24 <- fChange::ocidr(SPY_22_24)$data
colnames(data_22_24) <- colnames(SPY_22_24)

N_fit <- sum(
  as.Date(colnames(data_22_24), "%m-%d-%Y") < as.Date("2024-01-01")
)
results <- check_model_tightness(
  data = data_22_24,
  N_fit = N_fit
)

print(colMeans(results$viol_estimate))
print(colMeans(results$mean_estimate))

cols <- scales::hue_pal()(4)
png("./paper/figs/tightness1_2224.png", width = 800, height = 800)
print(
  ggplot2::ggplot() +
    ggplot2::geom_line(
      ggplot2::aes(
        x = seq_len(nrow(results$mean_estimate)),
        y = results$mean_estimate[, 2],
        color = "CCC-op-ARCH(1)"
      ), linewidth = 1
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        x = seq_len(nrow(results$mean_estimate)),
        y = results$mean_estimate[, 3],
        color = "CCC-op-ARCH(5)"
      ), linewidth = 1
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        x = seq_len(nrow(results$mean_estimate)),
        y = results$mean_estimate[, 4],
        color = "pw-fARCH(1)"
      ), linewidth = 1
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        x = seq_len(nrow(results$mean_estimate)),
        y = results$mean_estimate[, 5],
        color = "Historical"
      ), linewidth = 1
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = c(0.8, 0.9),
      axis.title = ggplot2::element_text(size = 30),
      axis.text = ggplot2::element_text(size = 24),
      legend.text = ggplot2::element_text(size = 24)
    ) +
    ggplot2::ylab("") +
    ggplot2::xlab("Resolution") +
    ggplot2::labs(color = NULL) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = max(1, floor(nrow(results$mean_estimate) * 0.5)),
        y = results$mean_estimate[
          max(1, floor(nrow(results$mean_estimate) * 0.5)),
          ] + c(0, 0.2, -0.1, 0.5, -0.3),
        label = c("", round(colMeans(results$viol_estimate), 3))
      ),
      size = 10
    ) +
    ggplot2::scale_color_manual(values = cols[c(1, 3, 4, 2)])
)
dev.off()

data("SPY_18_20", package = "fungarch")
data_18_20 <- fChange::ocidr(SPY_18_20)$data
colnames(data_18_20) <- colnames(SPY_18_20)
N_fit_18_20 <- sum(
  as.Date(colnames(data_18_20), "%m-%d-%Y") < as.Date("2020-01-01")
)
results_18_20 <- check_model_tightness(data_18_20, N_fit = N_fit_18_20)

print(colMeans(results_18_20$viol_estimate))
print(colMeans(results_18_20$mean_estimate))

png("./paper/figs/tightness1_1820.png", width = 800, height = 800)
print(
  ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(
      x = seq_len(nrow(results_18_20$mean_estimate)),
      y = results_18_20$mean_estimate[, 2], color = "CCC-op-ARCH(1)"
    ), linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(
      x = seq_len(nrow(results_18_20$mean_estimate)),
      y = results_18_20$mean_estimate[, 3], color = "CCC-op-ARCH(5)"
    ), linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(
      x = seq_len(nrow(results_18_20$mean_estimate)),
      y = results_18_20$mean_estimate[, 4], color = "pw-fARCH(1)"
    ), linewidth = 1) +
    ggplot2::geom_line(ggplot2::aes(
      x = seq_len(nrow(results_18_20$mean_estimate)),
      y = results_18_20$mean_estimate[, 5], color = "Historical"
    ), linewidth = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Resolution", y = "", color = NULL) +
    ggplot2::scale_color_manual(values = cols[c(1, 3, 4, 2)])
)
dev.off()
