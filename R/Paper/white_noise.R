library(devtools)
library(tidyverse)
##################################

SPY2018 <- read.csv("C:/Users/jerem/Downloads/SPY.USUSD_Candlestick_10_M_BID_01.01.2018-31.12.2018.csv")
SPY2019 <- read.csv("C:/Users/jerem/Downloads/SPY.USUSD_Candlestick_10_M_BID_01.01.2019-31.12.2019.csv")
SPY2020 <- read.csv("C:/Users/jerem/Downloads/SPY.USUSD_Candlestick_10_M_BID_01.01.2020-31.12.2020.csv")
# SPY2021 <- read.csv("C:/Users/jerem/Downloads/SPY.USUSD_Candlestick_10_M_BID_01.01.2021-31.12.2021.csv")
# SPY2022 <- read.csv("C:/Users/jerem/Downloads/SPY.USUSD_Candlestick_10_M_BID_01.01.2022-31.12.2022.csv")
# SPY2023 <- read.csv("C:/Users/jerem/Downloads/SPY.USUSD_Candlestick_10_M_BID_01.01.2023-31.12.2023.csv")
# SPY2024 <- read.csv("C:/Users/jerem/Downloads/SPY.USUSD_Candlestick_10_M_BID_01.01.2024-31.12.2024.csv")
# SPY2020 <- read.csv("C:/Users/jerem/Downloads/SPY.USUSD_Candlestick_10_M_BID_01.01.2020-31.12.2020.csv")

spy <- rbind(SPY2018, SPY2019, SPY2020)#,
             # SPY2021, SPY2022, SPY2023,
             # SPY2024)

spy$time <- as.POSIXct(spy$Local.time,
                       format = "%d.%m.%Y %H:%M:%OS ",
                       tz = "UTC")

#####
data = spy[,c("time","Close")]
# base drop
data <- data[as.numeric(format(data$time, '%H'))>=9 &
               as.numeric(format(data$time, '%H'))<=15, ]
data <- data[!(as.numeric(format(data$time, '%H'))==9 &
                 as.numeric(format(data$time, '%M'))<30),]
drop_vals = c()
dates = unique(format(data$time,'%Y-%m-%d'))
for(i in dates){
  cat(i,', ')
  if(length(unique(data[format(data$time,'%Y-%m-%d')==i,'Close']))==1){
    drop_vals = c(drop_vals, i)
  }
}
data <- data[!(format(data$time,'%Y-%m-%d') %in% drop_vals), ]
data$clock = format(data$time,'%H %M')
data$day = format(data$time,'%m-%d-%Y')
data1 <- data[,c('clock', 'day', 'Close')]
data_2224 <- data[, c('clock', 'day', 'Close')]
library(tidyverse)
plot_data <- data_2224 %>% pivot_wider(names_from = day,
                                       values_from = Close)
plot_data <- plot_data[,-1] %>% as.data.frame()

stock_wide = plot_data
for(i in 1:ncol(stock_wide)){
  # If any missing
  if(sum(is.na(stock_wide[,i]))>0){
    missing_idx <- which(is.na(stock_wide[,i]))

    nonmissing <- na.omit(stock_wide[,i])[[1]]
    first_val <- nonmissing[1]
    #lastval_val <- nonmissing[length(nonmissing)]

    for(j in missing_idx){
      if(j==1) stock_wide[j,i] <- first_val
      else stock_wide[j,i] <- stock_wide[j-1,i]
    }
  }
}

stock_wide <- as.data.frame(stock_wide)


#######################


errors <- function(resolution, alphas, basis, delta, data, C_eps, dimen){

  N <- ncol(data)

  #############
  # METHOD
  M <- array(dim=c(dimen,N))
  for(i in 1:N){
    for(j in 1:dimen){
      M[j,i] <- (data[,i] %*% basis[,j]) / nrow(data)
    }
  }


  # Get Alpha
  alpha_applied <- replicate(ncol(alphas),
                             matrix(0, nrow=resolution^2, ncol = resolution^2),
                             simplify=FALSE)
  for(p in 1:ncol(alphas)){
    for(k in 1:nrow(alphas)){
      hkk <- as.vector(basis[,k] %o% basis[,k])
      alphas_list <- sapply(1:ncol(alphas), function(d){
        alphas[k,d] * (hkk %o% hkk)
      },simplify = F)

      alpha_applied[[p]] <- alpha_applied[[p]] + Reduce('+', alphas_list)
    }
  }

  # Get Sigma
  sigma <- matrix(ncol=N, nrow=length(as.vector(delta)))
  for(i in (length(alpha_applied)+1):N){
    # sigma[,i] <- as.vector(delta) +
    #   operator(alpha_applied[[1]],
    #            tcrossprod(data[,i-1],data[,i-1])
    #   )

    sigma[,i] <- as.vector(delta) +
      Reduce('+',
             sapply(1:ncol(alphas),function(idx){
               operator(alpha_applied[[idx]],
                        tcrossprod(data[,i-idx],data[,i-idx])
               )
             },simplify = FALSE)
      )
  }


  errors <- matrix(nrow=dimen,ncol=N-length(alpha_applied))
  for(i in (length(alpha_applied)+1):N){
    sigmaP <- matrix(0, ncol=dimen, nrow=dimen)

    for(j in 1:dimen){
      for(k in 1:dimen){
        sigmaP[j,k] <-
          sum( matrix(sigma[,i],resolution,resolution) *
                 (basis[,j] %o% basis[,k]) ) / resolution^2
      }
    }

    # Get Sqrt(Sigma)
    sigma_svd <- svd(sigmaP)
    sigma_sqrt <-
      sigma_svd$u %*% diag(sqrt(sigma_svd$d)) %*% t(sigma_svd$v)

    errors[,i-length(alpha_applied)] <- solve(sigma_sqrt) %*% M[,i]
  }


  e_hat <- matrix(0, nrow=resolution, ncol=ncol(errors))
  for(i in 1:ncol(errors)){
    for(j in 1:dimen){
      e_hat[,i] <- e_hat[,i] +
        errors[j,i] * basis[,j]
    }
  }


  list('e_hat'=e_hat)
}

size_power <-
  function(data, w=0.05){

    resolution <- nrow(data)

    ##################
    ##  C_eps
    ##################
    ## Full C_eps
    C_eps <- fChange:::var(fChange:::dfts(data))#.generate_ceps( resolution )
    eig_Ce <- eigen( C_eps )
    basis <- eig_Ce$vectors * sqrt(resolution)
    eigenvalues <- eig_Ce$values / resolution


    ##################
    ##  Testing
    ##################
    X2 <- array(apply(data, 2, tcrossprod),
                dim=c(nrow(data),nrow(data),ncol(data)))

    ## Diag model with correct
    diag_mod <- .fit_diag_model_new(X = data,
                                    C_eps = C_eps,
                                    basis = basis,
                                    eigenvalues = eigenvalues,
                                    w=w, w1=w, p=1)
    # PREDICTIONS
    errors_tmp = errors(resolution, alphas = diag_mod$alpha,
                        basis=basis, delta=diag_mod$delta,
                        data=data, C_eps=C_eps, dimen=diag_mod$K)

    ## Diag 5
    diag_mod5 <- .fit_diag_model_new(X = data,
                                     C_eps = C_eps,
                                     basis = basis,
                                     eigenvalues = eigenvalues,
                                     w=w, w1=w, p=5)
    # PREDICTIONS
    errors5_tmp = errors(resolution, alphas = diag_mod5$alpha,
                        basis=basis, delta=diag_mod5$delta,
                        data=data, C_eps=C_eps, dimen=diag_mod5$K)


    ########
    N <- ncol(data)
    M <- array(dim=c(resolution,N))
    for(i in 1:N){
      for(j in 1:resolution){
        M[j,i] <- (data[,i] %*% basis[,j]) / nrow(data)
      }
    }

    x_hat <- matrix(0, nrow=resolution, ncol=N)
    for(i in 1:N){
      for(j in 1:resolution){
        x_hat[,i] <- x_hat[,i] +
          M[j,i] * basis[,j]
      }
    }

    list('x_hat'=x_hat,
         'e_hat'=errors_tmp$e_hat,
         'e_hat5'=errors5_tmp$e_hat)
  }

load_all()

##################
data <- compute_ocidr(stock_wide)
# saveRDS(stock_wide, "C:/Users/jerem/Downloads/stock_wide.rds")
# saveRDS(data, "C:/Users/jerem/Downloads/stock_wide_data.rds")
stock_wide <- readRDS("C:/Users/jerem/Downloads/stock_wide.rds")
data <- readRDS("C:/Users/jerem/Downloads/stock_wide_data.rds")

set.seed(1234)
res <- size_power(data = data)


fChange::portmanteau_tests(res$x_hat^2,test = 'multi',lag=3)$pvalue
fChange::portmanteau_tests(res$e_hat^2,test = 'multi',lag=3)$pvalue
fChange::portmanteau_tests(res$e_hat5^2,test = 'multi', lag=3)$pvalue

fChange::portmanteau_tests(res$x_hat^2,test = 'multi',lag=10)$pvalue
fChange::portmanteau_tests(res$e_hat^2,test = 'multi',lag=10)$pvalue
fChange::portmanteau_tests(res$e_hat5^2,test = 'multi', lag=10)$pvalue


#################
# install.packages("reticulate")
Sys.setenv(RETICULATE_PYTHON = "C:/Users/jerem/AppData/Local/Programs/Python/Python314/python.exe")
library(reticulate)
library(tidyverse)
# use_python(install_python())
# py_install(c("kaleido", "plotly"))
# py_install('kaleido')

# library(reticulate)
# reticulate::use_python("C:/Users/jerem/AppData/Local/Programs/Python/Python314/python.exe",
#            required = TRUE)
# # reticulate::install_python()
# reticulate::py_run_string("import sys")
# #pip install --upgrade "kaleido==0.1.*"

# saveRDS(res$e_hat, 'C:/Users/jerem/Downloads/tmp1.rds')
# saveRDS(as.Date(colnames(stock_wide),"%m-%d-%Y"),'C:/Users/jerem/Downloads/names.rds')
X <- fChange::dfts(res$e_hat,
                   labels = as.Date(colnames(stock_wide),"%m-%d-%Y")[-1],
                   inc.warnings = FALSE)

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

# devtools::install_github('quinten-goens/plotly.R@fix/kaleido-export-bug')
plotly::save_image(
  plt,
  file = "C:/Users/jerem/Downloads/sp500_resids_1_full.png",
)
plotly::save_image(
  plt,
  file = "C:\\Users\\jerem\\Downloads\\sp500_resids_1_full.png",
)

#################

# saveRDS(res$e_hat, 'C:/Users/jerem/Downloads/tmp1.rds')
# saveRDS(as.Date(colnames(stock_wide),"%m-%d-%Y"),'C:/Users/jerem/Downloads/names.rds')
X <- fChange::dfts(res$e_hat5,
                   labels = as.Date(colnames(stock_wide),"%m-%d-%Y")[-c(1:5)],
                   inc.warnings = FALSE)

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

plotly::save_image(
  plt,
  file = "C:/Users/jerem/Downloads/sp500_resids_5.png"
)

##################


png('C:/Users/jerem/Downloads/x_orig_wn.png',width = 800, height =800)
fChange::sacf(data, lag.max = 10)$plot +
  theme(axis.text = element_text(size=40),
        axis.title = element_text(size=50))
dev.off()
png('C:/Users/jerem/Downloads/x_hat_orig_wn.png',width = 800, height =800)
fChange::sacf(res$x_hat, lag.max = 10)$plot +
  theme(axis.text = element_text(size=40),
        axis.title = element_text(size=50))
dev.off()

png('C:/Users/jerem/Downloads/x_hat_wn.png',width = 800, height =800)
fChange::sacf(res$x_hat^2, lag.max = 10)$plot +
  theme(axis.text = element_text(size=40),
        axis.title = element_text(size=50))
dev.off()

png('C:/Users/jerem/Downloads/e_hat_wn.png',width = 800, height =800)
fChange::sacf(res$e_hat^2, lag.max = 10)$plot +
  theme(axis.text = element_text(size=40),
        axis.title = element_text(size=50))
dev.off()

png('C:/Users/jerem/Downloads/e_hat5_wn.png',width = 800, height =800)
fChange::sacf(res$e_hat5^2, lag.max = 10)$plot +
  theme(axis.text = element_text(size=40),
        axis.title = element_text(size=50))
dev.off()






