
check_onestep <- function(data, start = 500, est_dimen=3, w=0.1, silent=FALSE){
  ##################
  ## Params
  ##################
  resolution <- nrow(data)
  N <- ncol(data)

  C_eps <- .generate_ceps(resolution)
  eig_Ce <- eigen( C_eps )
  basis <- eig_Ce$vectors * sqrt(resolution)
  eigenvalues <- eig_Ce$values / resolution

  ##################
  ## CV
  ##################
  X2 <- array(apply(data, 2, tcrossprod),
              dim=c(resolution, resolution,N))
  errors <- matrix(NA, ncol=9,nrow=N-start)
  colnames(errors) <- c('diag','diag2','pw', 'pw-t',
                        'scalar','scalar1',
                        'last','mean','zero')
  for(i in (start+1):N){
    if(!silent) cat(i-start,', ')
    X2_prev <- X2[,,i-1]

    ## Diag model with correct
    diag_mod <- .fit_diag_model(X = data[,1:(i-1)],
                                C_eps = C_eps,
                                basis = basis,
                                eigenvalues = eigenvalues,
                                w=w, p=1, dimen=est_dimen)
    # TODO:: Delta Information
    diag_pred <- .pred_diag(resolution, alphas = diag_mod$alpha,
                            basis=basis, delta=diag_mod$delta, X2_prev=X2_prev,
                            C_eps=C_eps)
    # diag_pred1 <- .pred_diag(resolution, alphas = diag_mod$alpha,
    #                          basis=basis, delta=diag_mod$delta1, X2_prev=X2_prev,
    #                          C_eps=C_eps)
    diag_pred2 <- .pred_diag(resolution, alphas = diag_mod$alpha,
                             basis=basis, delta=diag_mod$delta2, X2_prev=X2_prev,
                             C_eps=C_eps)


    ## PW Model - MP
    pw_mod <- .fit_pw(X=data[,1:(i-1)], w=w)
    pw_pred <- .pred_pw(X=data[,1:(i-1)], pw_mod$alpha,
                        pw_mod$delta, C_eps)

    ## PW Model - Tik
    pw_mod1 <- pointwiseARCH(X=data[,1:(i-1)],w = 0.001)
    pw_pred1 <- .pred_pw(X=data[,1:(i-1)], pw_mod1$hat_alpha,
                         pw_mod1$hat_delta, C_eps)

    ## scalar
    sink(nullfile())
    arch1_pred <- tryCatch({
      arch1_mod <- fGarch::garchFit(formula = ~ arma(0,0) + garch(1,0),
                                    data = as.numeric(data[1,])
      )
      fGarch::predict(arch1_mod)$standardDeviation[1] *
        matrix(1,nrow(X2),nrow(X2))
    }, error=function(e){
      matrix(NA,nrow(X2),nrow(X2))
    })

    arch1_mod1 <- fGarch::garchFit(formula = ~ arma(0,0) + garch(1,0),
                                  data = as.numeric(colMeans(data))
    )
    arch1_pred1 <- fGarch::predict(arch1_mod1)$standardDeviation[1] *
      matrix(1,nrow(X2),nrow(X2))
    sink()

    ## Last Value
    last_pred <- X2[,,i-1]
    ## Mean
    mean_pred <- rowMeans(X2[,,1:(i-1)], dims = 2)
    ## Zero
    zero_pred <- matrix(0,resolution,resolution)


    # Compute Errors
    #   TODO:: Verify this
    errors[i-start,] <- c(
      dot_integrate( (X2[,,i] - diag_pred/N)^2 ),
      #dot_integrate( (X2[,,i] - diag_pred1/N)^2 ),
      dot_integrate( (X2[,,i] - diag_pred2/N)^2 ),
      dot_integrate( (X2[,,i] - pw_pred)^2 ),
      dot_integrate( (X2[,,i] - pw_pred1)^2 ),
      dot_integrate( (X2[,,i] - arch1_pred)^2 ),
      dot_integrate( (X2[,,i] - arch1_pred1)^2 ),
      dot_integrate( (X2[,,i] - last_pred)^2 ),
      dot_integrate( (X2[,,i] - mean_pred)^2 ),
      dot_integrate( (X2[,,i] - zero_pred)^2 )
    )
  }

  errors
}

check_onestep1 <- function(data, start = 500, est_dimen=3, w=0.1, silent=FALSE){
  ##################
  ## Params
  ##################
  resolution <- nrow(data)
  N <- ncol(data)

  C_eps <- .generate_ceps(resolution)
  eig_Ce <- eigen( C_eps )
  basis <- eig_Ce$vectors * sqrt(resolution)
  eigenvalues <- eig_Ce$values / resolution

  ##################
  ## CV
  ##################
  X2 <- array(apply(data, 2, tcrossprod),
              dim=c(resolution, resolution,N))
  errors <- matrix(NA, ncol=9,nrow=N-start)
  colnames(errors) <- c('diag','diag2','pw', 'pw-t',
                        'scalar','scalar1',
                        'last','mean','zero')
  for(i in (start+1):N){
    if(!silent) cat(i-start,', ')
    X2_prev <- X2[,,i-1]

    ## Diag model with correct
    diag_mod <- .fit_diag_model(X = data[,1:(i-1)],
                                C_eps = C_eps,
                                basis = basis,
                                eigenvalues = eigenvalues,
                                w=w, p=1, dimen=est_dimen)
    # TODO:: Delta Information
    diag_pred <- .pred_diag(resolution, alphas = diag_mod$alpha,
                            basis=basis, delta=diag_mod$delta, X2_prev=X2_prev,
                            C_eps=C_eps)
    # diag_pred1 <- .pred_diag(resolution, alphas = diag_mod$alpha,
    #                          basis=basis, delta=diag_mod$delta1, X2_prev=X2_prev,
    #                          C_eps=C_eps)
    diag_pred2 <- .pred_diag(resolution, alphas = diag_mod$alpha,
                             basis=basis, delta=diag_mod$delta2, X2_prev=X2_prev,
                             C_eps=C_eps)


    ## PW Model - MP
    pw_mod <- .fit_pw(X=data[,1:(i-1)], w=w)
    pw_pred <- .pred_pw(X=data[,1:(i-1)], pw_mod$alpha,
                        pw_mod$delta, C_eps)

    ## PW Model - Tik
    pw_mod1 <- pointwiseARCH(X=data[,1:(i-1)],w = 0.001)
    pw_pred1 <- .pred_pw(X=data[,1:(i-1)], pw_mod1$hat_alpha,
                         pw_mod1$hat_delta, C_eps)

    ## scalar
    sink(nullfile())
    arch1_pred <- tryCatch({
      arch1_mod <- fGarch::garchFit(formula = ~ arma(0,0) + garch(1,0),
                                    data = as.numeric(data[1,])
      )
      fGarch::predict(arch1_mod)$standardDeviation[1] *
        matrix(1,nrow(X2),nrow(X2))
    }, error=function(e){
      matrix(NA,nrow(X2),nrow(X2))
    })

    arch1_mod1 <- fGarch::garchFit(formula = ~ arma(0,0) + garch(1,0),
                                   data = as.numeric(colMeans(data))
    )
    arch1_pred1 <- fGarch::predict(arch1_mod1)$standardDeviation[1] *
      matrix(1,nrow(X2),nrow(X2))
    sink()

    ## Last Value
    last_pred <- X2[,,i-1]
    ## Mean
    mean_pred <- rowMeans(X2[,,1:(i-1)], dims = 2)
    ## Zero
    zero_pred <- matrix(0,resolution,resolution)


    # Compute Errors
    #   TODO:: Verify this
    errors[i-start,] <- c(
      dot_integrate( abs(X2[,,i] - diag_pred/N) ),
      #dot_integrate( abs(X2[,,i] - diag_pred1/N) ),
      dot_integrate( abs(X2[,,i] - diag_pred2/N) ),
      dot_integrate( abs(X2[,,i] - pw_pred) ),
      dot_integrate( abs(X2[,,i] - pw_pred1) ),
      dot_integrate( abs(X2[,,i] - arch1_pred) ),
      dot_integrate( abs(X2[,,i] - arch1_pred1) ),
      dot_integrate( abs(X2[,,i] - last_pred) ),
      dot_integrate( abs(X2[,,i] - mean_pred) ),
      dot_integrate( abs(X2[,,i] - zero_pred) )
    )
  }

  errors
}


compute_cidr <- function(dat){
  dat_cidr <- dat
  for(i in 1:nrow(dat_cidr)){
    dat_cidr[i,] <- 100*(log(dat[i,]) - log(dat[1,]))
  }
  dat_cidr
}

compute_ocidr <- function(dat){
  dat_cidr <- dat
  dat_cidr[,1]  <- 100*(log(dat[,1]) - log(dat[1,1]))
  # Obs
  for(j in 2:ncol(dat_cidr)){
    dat_cidr[,j] <- 100 *
      ( log(dat[,j]) - log(dat[nrow(dat_cidr),j-1]) )

    # for(i in 2:nrow(dat_cidr)){
    #   dat_cidr[i,j] <- 100*
    #     ( log(dat[i,j]) - log(dat[1,j]) )
    # }
  }
  dat_cidr
}

compute_onr <- function(dat){
  dat_nr <- dat
  dat_nr[1,1] <- 0
  dat_nr[1,-1]  <- 100 * ( log(dat[1,-1]) - log(dat[1,-ncol(dat)]) )
  # Obs
  for(j in 1:ncol(dat_nr)){
    dat_nr[-1,j] <- 100 *
      ( log(dat[-1,j]) - log(dat[-nrow(dat_nr),j]) )
  }
  dat_nr
}

compute_nr <- function(dat){
  dat_nr <- dat
  dat_nr[1,] <- 0
  # Obs
  for(j in 1:ncol(dat_nr)){
    dat_nr[-1,j] <- 100 *
      ( log(dat[-1,j]) - log(dat[-nrow(dat_nr),j]) )
  }
  dat_nr
}

getData <- function(company){
  year2021 <-
    read.csv(paste0("C:/Users/jerem/Downloads/stocks/",
                    company,"2021.csv"))
  year2022 <-
    read.csv(paste0("C:/Users/jerem/Downloads/stocks/",
                    company,"2022.csv"))
  year2023 <-
    read.csv(paste0("C:/Users/jerem/Downloads/stocks/",
                    company,"2023.csv"))

  stock <- rbind(year2021,
                 year2022,
                 year2023)
  stock <- stock[stock$Volume!=0,]
  stock$datetime <-
    as.POSIXct(stock$Local.time,
               format="%d.%m.%Y %H:%M:%S")
  stock$date <- as.Date(stock$datetime)
  stock$hourMin <- paste(hour(stock$datetime),'-',minute(stock$datetime))
  stock$hr <- hour(stock$datetime)
  stock$min <- minute(stock$datetime)
  stock_wide <- stock %>%
    pivot_wider(names_from = date,
                values_from = Close,
                id_cols = c(hr,min) ) %>%
    arrange(.,hr,min)

  for(i in 3:ncol(stock_wide)){
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

  stock_wide <- as.data.frame(stock_wide[,-c(1:2)])

  stock_wide %>%
    mutate(grp = 1+ (row_number()-1) %/% 15) %>%
    group_by(grp) %>%
    # summarise(across(everything(), mean, na.rm = TRUE)) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%
    select(-grp) %>%
    as.data.frame()
}


#################

library(devtools)
library(tidyverse)
load_all()

msft_15min <- getData("MSFT")
msft_cidr15 <- compute_cidr(msft_15min)
msft_ocidr15 <- compute_ocidr(msft_15min)
msft_nr15 <- compute_nr(msft_15min)
msft_onr15 <- compute_onr(msft_15min)

res_msft0 <- check_onestep(msft_cidr15)
res_msft1 <- check_onestep(msft_ocidr15)
res_msft2 <- check_onestep(msft_nr15)
res_msft3 <- check_onestep(msft_onr15)
res_msft0a <- check_onestep1(msft_cidr15)
res_msft1a <- check_onestep1(msft_ocidr15)
res_msft2a <- check_onestep1(msft_nr15)
res_msft3a <- check_onestep1(msft_onr15)

amzn_15min <- getData("AMZN")
amzn_cidr15 <- compute_cidr(amzn_15min)
amzn_ocidr15 <- compute_ocidr(amzn_15min)
amzn_nr15 <- compute_nr(amzn_15min)
amzn_onr15 <- compute_onr(amzn_15min)

res_amzn0 <- check_onestep(amzn_cidr15)
res_amzn1 <- check_onestep(amzn_ocidr15)
res_amzn2 <- check_onestep(amzn_nr15)
res_amzn3 <- check_onestep(amzn_onr15)
res_amzn0a <- check_onestep1(amzn_cidr15)
res_amzn1a <- check_onestep1(amzn_ocidr15)
res_amzn2a <- check_onestep1(amzn_nr15)
res_amzn3a <- check_onestep1(amzn_onr15)

tsla_15min <- getData("TSLA")
tsla_cidr15 <- compute_cidr(tsla_15min)
tsla_ocidr15 <- compute_ocidr(tsla_15min)
tsla_nr15 <- compute_nr(tsla_15min)
tsla_onr15 <- compute_onr(tsla_15min)

res_tsla0 <- check_onestep(tsla_cidr15)
res_tsla1 <- check_onestep(tsla_ocidr15)
res_tsla2 <- check_onestep(tsla_nr15)
res_tsla3 <- check_onestep(tsla_onr15)
res_tsla0a <- check_onestep1(tsla_cidr15)
res_tsla1a <- check_onestep1(tsla_ocidr15)
res_tsla2a <- check_onestep1(tsla_nr15)
res_tsla3a <- check_onestep1(tsla_onr15)

colMeans(res_msft0)
colMeans(res_amzn0)
colMeans(res_tsla0)
colMeans(res_msft1)
colMeans(res_amzn1)
colMeans(res_tsla1)
colMeans(res_msft2)
colMeans(res_amzn2)
colMeans(res_tsla2)
colMeans(res_msft3)
colMeans(res_amzn3)
colMeans(res_tsla3)

colMeans(res_msft0a)
colMeans(res_amzn0a)
colMeans(res_tsla0a)
colMeans(res_msft1a)
colMeans(res_amzn1a)
colMeans(res_tsla1a)
colMeans(res_msft2a)
colMeans(res_amzn2a)
colMeans(res_tsla2a)
colMeans(res_msft3a)
colMeans(res_amzn3a)
colMeans(res_tsla3a)

# window/re-estimate
# abs
# cidr
