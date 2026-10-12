#####
# Both Examples
N <- 200
resolution <- 50

#####
# Example 1

## C_eps
C_eps <- .generate_ceps(resolution,denom = 1)

eig_Ce <- eigen( C_eps )
basis <- eig_Ce$vectors * sqrt(resolution)
eigenvalues <- eig_Ce$values / resolution
C_eps_fit <- matrix(0, nrow=resolution, ncol = resolution)

alpha_scale <- c(0,rep(1.6,3))/eigenvalues[1:4]

## SETUP ALPHA/DELTA
#   Alpha is based on tensor sum of `alpha_dimen` basis functions
alpha <- matrix(0, nrow=resolution^2, ncol = resolution^2)
for(k in 1:length(alpha_scale)){
  hkk <- as.vector(basis[,k] %o% basis[,k])
  alpha <- alpha + alpha_scale[k] * (hkk %o% hkk)
}

#   Delta is full C_eps
delta <- .generate_ceps(resolution)


## Generate Data
set.seed(123)
data <- generateFARCH_use(N = N, resolution = resolution,
                          alpha=alpha, delta=delta,
                          basis=basis,
                          burnin=100, p=1,
                          eigs=eigenvalues[1:length(alpha_scale)],#eigs_gen,
                          C_eps=C_eps)

X <- fChange::dfts(data$data)
X <- fChange::dfts(X, inc.warnings = FALSE)

changes = NULL
plot_title = NULL
val_axis_title = ""
res_axis_title = ""
FD_axis_title = ""
eye = list(x = 0.25, y = -1.5, z = 0.1)
aspectratio = list(x = 1.1, y = 0.75, z = 0.75)
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
                 showticklabels = FALSE,
                 ticktext = fChange:::.select_n(vals = X$labels, n = 6),
                 tickvals = fChange:::.select_n(vals = 1:length(X$labels), n = 6)),
    zaxis = list(title = val_axis_title,
                 tickfont = list(size = 15),
                 showticklabels = showticklabels))) %>%
  plotly::layout(title = plot_title, scene = scene) %>%
  plotly::layout(showlegend = FALSE)

plotly::save_image(
  plt,
  file = "./paper/example1.png"
)

#####
# Example 2

##  C_eps
C_eps <- matrix(NA, resolution, resolution)
v <- seq(0, 1, length.out = resolution)
for (i in 1:resolution) {
  for(j in 1:resolution){
    C_eps[i,j] <- min(v[j],v[i])
  }
}

eig_Ce <- eigen( C_eps )
basis <- eig_Ce$vectors * sqrt(resolution)
eigenvalues <- eig_Ce$values / resolution
C_eps_fit <- matrix(0, nrow=resolution, ncol = resolution)

alpha_scale <- c(0,rep(1.1,3))/eigenvalues[1:4]

## SETUP ALPHA/DELTA
#   Alpha is based on tensor sum of `alpha_dimen` basis functions
alpha <- matrix(0, nrow=resolution^2, ncol = resolution^2)
for(k in 1:length(alpha_scale)){
  hkk <- as.vector(basis[,k] %o% basis[,k])
  alpha <- alpha + alpha_scale[k] * (hkk %o% hkk)
}

#   Delta is full C_eps
delta <- C_eps

## Generate Data
set.seed(1234567)
data <- generateFARCH_use(N = N, resolution = resolution,
                          alpha=alpha,
                          delta=delta +
                            rnorm(resolution,sd=1/7) * delta,#,rep(0,resolution-4))),
                          basis=basis,
                          burnin=100, p=1,
                          eigs=eigenvalues[1:length(alpha_scale)],#eigs_gen,
                          C_eps=C_eps)

X <- fChange::dfts(data$data, inc.warnings = FALSE)

changes = NULL
plot_title = NULL
val_axis_title = ""
res_axis_title = ""
FD_axis_title = ""
eye = list(x = 0.25, y = -1.5, z = 0.1)
aspectratio = list(x = 1.1, y = 0.75, z = 0.75)
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
                 showticklabels = FALSE,
                 ticktext = fChange:::.select_n(vals = X$labels, n = 6),
                 tickvals = fChange:::.select_n(vals = 1:length(X$labels), n = 6)),
    zaxis = list(title = val_axis_title,
                 tickfont = list(size = 15),
                 showticklabels = showticklabels))) %>%
  plotly::layout(title = plot_title, scene = scene) %>%
  plotly::layout(showlegend = FALSE)

plotly::save_image(
  plt,
  file = "./paper/example2.png"
)
