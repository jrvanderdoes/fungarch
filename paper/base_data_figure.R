#####
# Get Data
X <- fChange::dfts(spy,labels = as.Date(colnames(spy),"%m-%d-%Y"))

#####
# Base Plot
changes = NULL
plot_title = NULL
val_axis_title = ""
res_axis_title = ""
FD_axis_title = ""
eye = list(x = -0.5, y = -1.5, z = 0.25)
aspectratio = list(x = 1.1, y = 0.5, z = 0.5)
showticklabels = TRUE

## PLOTTING FUNCTION
plotData <- data.frame()
for (i in 1:ncol(X)) {
  plotData <- rbind(plotData, data.frame(resolution = X$fparam,
                                         FDRep = i, Value = X$data[, i]))
}
scene <- list(camera = list(eye = eye),
              aspectmode = "manual",
              aspectratio = aspectratio)
plot_colors <- RColorBrewer::brewer.pal(11, "Spectral")
plot_colors[6] <- "yellow"
plt <- plotly::plot_ly(plotData, x = ~FDRep,
                       y = ~resolution,
                       z = ~Value,
                       type = "scatter3d",
                       mode = "lines",
                       color = ~as.factor(FDRep),
                       colors = plot_colors) %>%
  plotly::layout(scene = list(
    yaxis = list(title = res_axis_title,
                 showticklabels = FALSE,
                 ticktext = round(seq(min(X$fparam),
                                      max(X$fparam), length.out = 5), 3),
                 tickvals = fChange:::.select_n(vals = seq(min(X$fparam),
                                                           max(X$fparam), length.out = 5), n = 5)),
    xaxis = list(title = FD_axis_title,
                 tickfont = list(size = 15),
                 showticklabels = showticklabels,
                 ticktext = fChange:::.select_n(vals = format(as.Date(X$labels), "%b %y"), n = 6),
                 tickvals = fChange:::.select_n(vals = 1:length(X$labels), n = 6)),
    zaxis = list(title = val_axis_title,
                 showticklabels = FALSE))) %>%
  plotly::layout(title = plot_title, scene = scene) %>%
  plotly::layout(showlegend = FALSE)

#####
# Export Base Plot
plotly::save_image(
  plt,
  file = "./paper/sp500.png"
)


#####
# Get Data
# X <- fChange::dfts(compute_ocidr(spy),
#                    labels = as.Date(colnames(spy),"%m-%d-%Y"))
X <- fChange::ocidr(spy)
X$labels <- as.Date(colnames(spy),"%m-%d-%Y")

#####
# OCIDR Plot
changes = NULL
plot_title = NULL
val_axis_title = ""
res_axis_title = ""
FD_axis_title = ""
eye = list(x = -0.5, y = -1.5, z = 0.25)
aspectratio = list(x = 1.1, y = 0.5, z = 0.5)
showticklabels = TRUE

## PLOTTING FUNCTION
plotData <- data.frame()
for (i in 1:ncol(X)) {
  plotData <- rbind(plotData, data.frame(resolution = X$fparam,
                                         FDRep = i, Value = X$data[, i]))
}
scene <- list(camera = list(eye = eye),
              aspectmode = "manual",
              aspectratio = aspectratio)
plot_colors <- RColorBrewer::brewer.pal(11, "Spectral")
plot_colors[6] <- "yellow"
plt <- plotly::plot_ly(plotData, x = ~FDRep,
                       y = ~resolution,
                       z = ~Value,
                       type = "scatter3d",
                       mode = "lines",
                       color = ~as.factor(FDRep),
                       colors = plot_colors) %>%
  plotly::layout(scene = list(
    yaxis = list(title = res_axis_title,
                 showticklabels = FALSE,
                 ticktext = round(seq(min(X$fparam),
                                      max(X$fparam), length.out = 5), 3),
                 tickvals = fChange:::.select_n(vals = seq(min(X$fparam),
                                                           max(X$fparam), length.out = 5), n = 5)),
    xaxis = list(title = FD_axis_title,
                 tickfont = list(size = 15),
                 showticklabels = showticklabels,
                 ticktext = fChange:::.select_n(vals = format(as.Date(X$labels), "%b %y"), n = 6),
                 tickvals = fChange:::.select_n(vals = 1:length(X$labels), n = 6)),
    zaxis = list(title = val_axis_title,
                 showticklabels = FALSE))) %>%
  plotly::layout(title = plot_title, scene = scene) %>%
  plotly::layout(showlegend = FALSE)


#####
# Export OCIDR Plot
plotly::save_image(
  plt,
  file = "./paper/sp500_ocidr.png"
)
