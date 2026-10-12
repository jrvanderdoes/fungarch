devtools::load_all()
data("SPY_18_20", package = "fungarch")
dir.create("./paper/figs", recursive = TRUE, showWarnings = FALSE)

dat <- fChange::ocidr(SPY_18_20)$data
colnames(dat) <- colnames(SPY_18_20)

fit_amt <- sum(as.Date(colnames(dat), "%m-%d-%Y") < as.Date("2020-02-20"))
est_amt <- sum(as.Date(colnames(dat), "%m-%d-%Y") < as.Date("2020-03-21")) - fit_amt

res1 <- forecast_envelope(
  data = dat,
  N_fit = fit_amt,
  forecast = est_amt
)

png("./paper/figs/forecast_ocidr_spy1820.png",
    width = 1200, height = 400)
res1$plt
dev.off()

python_path <- Sys.getenv("RETICULATE_PYTHON")
if (nzchar(python_path)) {
  reticulate::use_python(python_path, required = TRUE)
}
# pip install --upgrade "kaleido==0.1.*"
# ? pip install plotly==5.3.1
# Also try pip install kaleido==0.2.1

plt_save <- plot_forecast_envelope(
  fChange:::dfts(dat[, 1:(fit_amt + est_amt)]),
  lower = res1$lower,
  upper = res1$upper,
  plot_title = NULL,
  eye = list(x = 0.75, y = -1.5, z = 0.75),
  val_axis_title = "", res_axis_title = "", FD_axis_title = "",
  showticklabels = FALSE
)
plotly::save_image(
  plt_save,
  file = "./paper/figs/spy_forecast_ocidr_spy1820.png"
)

plt_save <- plot_forecast_envelope(
  fChange:::dfts(dat[, fit_amt:(fit_amt + est_amt)]),
  lower = res1$lower,
  upper = res1$upper,
  plot_title = NULL,
  eye = list(x = 2, y = 0, z = 0.1),
  aspectratio = list(x = 2, y = 0.5, z = 0.5),
  val_axis_title = "", res_axis_title = "", FD_axis_title = "",
  showticklabels = FALSE, fixme = FALSE
)
plotly::save_image(
  plt_save,
  file = "./paper/figs/spy_forecast_close_ocidr_spy1820.png"
)
