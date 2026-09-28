#' Title
#'
#' @param resolution
#' @param N_fit
#' @param N_pred
#' @param alpha
#' @param delta
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
generateFitCompareOPARCH <-
  function(
    resolution = 24, N_fit = 500, N_pred = 250,
    alpha = function(s,t,u,v){ s*(1-s)*t*(1-t)*u*(1-u)*v*(1-v) },
    delta = function(s,t){ s*(1-s)*t*(1-t) },
    seed=NA){
    if(!is.na(seed)) set.seed(seed)

    dat <- generateFunctionalCondVol_ff(resolution = resolution,
                                        N = N_fit+N_pred,'arch',
                                        alpha_par = alpha,
                                        delta_par = delta)
    dat_fit <- dat$garch_mat[,1:N_fit]
    dat_pred <- dat$garch_mat[,N_fit+ (0:N_pred)]

    ######
    ## MLE
    ######
    NBASIS = 20
    fd = fda::Data2fd(argvals=seq(0,1,len=nrow(dat$garch_mat)),
                      y=dat_fit,
                      fda::create.bspline.basis(nbasis=NBASIS))

    # extract data-driven basis functions through the truncated FPCA method.
    basis_est = basis.est(dat_fit, M=2, "tfpca",NBASIS=NBASIS)$basis

    # estimate an FARCH(1) model with basis when M=1.
    arch1_est = suppressMessages( est.fArch(fd, basis_est[,1]) )
    alpha_est <- as.numeric(arch1_est$As[[1]]) * (as.matrix(basis_est[,1]) %*% t(as.matrix(basis_est[,1])))
    delta_est <- basis_est[,1] * arch1_est$d

    result_mle <- list('hat_alpha'=alpha_est,
                       'hat_delta'=delta_est)

    ######
    ## Pointwise
    ######
    result_pw <- pointwiseARCH(X = dat_fit, w = 0.01)

    ######
    ## Operator - YW
    ######
    result_yw <- operatorARCH(X = dat_fit, w = 0.01,p = 1 )

    ######
    ## Operator - MP
    ######
    result_mp90 <- operatorARCHmp(X = dat_fit, w = 0.1,p = 1 )
    result_mp95 <- operatorARCHmp(X = dat_fit, w = 0.05,p = 1 )
    result_mp99 <- operatorARCHmp(X = dat_fit, w = 0.01,p = 1 )

    ######
    ## Get error for all
    ######

    # C_eps
    C_eps <- matrix(NA,resolution,resolution)
    times <- 1:resolution/resolution
    for (i in 1:resolution){
      #C_eps[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
      C_eps[i,] <- exp(-abs(times[i]-times)/2)
    }

    dat_pred2_pw <- dat_pred^2
    dat_pred2 <- array(apply(dat_pred, 2, tcrossprod),
                       dim=c(nrow(dat_pred),nrow(dat_pred),ncol(dat_pred)))
    estimate_or <-
      .predict_ARCH_op_or(X2 = dat_pred2,
                          sigma_sqrt_list=dat$params$Sigma2_sqrt_list[N_fit+0:N_pred],
                          C_eps = C_eps, p = 1)

    estimate_yw <-
      .predict_ARCH_op(X2=dat_pred2,
                       alpha=result_yw$hat_alpha,
                       delta=result_yw$hat_delta,
                       C_eps=C_eps, p=1)

    estimate_mp90 <-
      .predict_ARCH_op(X2=dat_pred2,
                       alpha=result_mp90$hat_alpha,
                       delta=result_mp90$hat_delta,
                       C_eps=C_eps, p=1)
    estimate_mp95 <-
      .predict_ARCH_op(X2=dat_pred2,
                       alpha=result_mp95$hat_alpha,
                       delta=result_mp95$hat_delta,
                       C_eps=C_eps, p=1)
    estimate_mp99 <-
      .predict_ARCH_op(X2=dat_pred2,
                       alpha=result_mp99$hat_alpha,
                       delta=result_mp99$hat_delta,
                       C_eps=C_eps, p=1)

    estimate_pw <-
      .pointwise_next_step_ests_op(dat_pred2_pw,
                                   result_pw$hat_alpha,
                                   result_pw$hat_delta,
                                   C_eps)
    estimate_mle <-
      .pointwise_next_step_ests_op(dat_pred2_pw,
                                   result_mle$hat_alpha,
                                   result_mle$hat_delta,
                                   C_eps)

    ## Get Errors
    full_dat2 <- array(
      apply(dat$garch_mat, 2, tcrossprod),
      dim=c(nrow(dat$garch_mat),nrow(dat$garch_mat),ncol(dat$garch_mat)))
    error <-
      data.frame('OR'=rep(0, ncol(dat_pred2)),
                 'YW'=0,
                 'MP90'=0,'MP95'=0,'MP99'=0,
                 'PW'=0,
                 'MLE'=0,
                 'Mean'=0,
                 'Prev'=0,
                 'Zero'=0)
    for(i in 2:dim(dat_pred2)[3]){

      # Errors for next step (learning upto last, no re-fit)
      error[i,] <- c(
        dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_or[,,i])^2)),
        dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_yw[,,i])^2)),
        dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_mp90[,,i])^2)),
        dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_mp95[,,i])^2)),
        dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_mp99[,,i])^2)),
        dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_pw[,,i])^2)),
        dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_mle[,,i])^2)),
        dot_integrate(dot_integrate_col((dat_pred2[,,i] -
                                           rowMeans(full_dat2[,,1:(N_fit+i-2)],dims=2))^2)),
        dot_integrate(dot_integrate_col((dat_pred2[,,i] - dat_pred2[,,i-1])^2)),
        dot_integrate(dot_integrate_col((dat_pred2[,,i] - 0)^2))
      )
    }

    colSums(error)
  }

.pointwise_next_step_ests_op <- function(X2, alpha, delta, C_eps){

  estimate <- .predictARCH(X2, alpha, delta, est2=T)


  Estimate <- array(NA, dim=c(nrow(X2),nrow(X2),ncol(X2)))
  for(i in 1:ncol(estimate)){
    sqrt_est <- sqrt(pmax(estimate[,i],0))
    Estimate[,,i] <- t(sqrt_est %*% C_eps) %*% sqrt_est / length(sqrt_est)
  }

  Estimate
}

.predict_ARCH_op_or <- function(X2, sigma_sqrt_list, C_eps, p){
  # dims <- dim(X2)

  # Eq 3 (Eq 29)
  # Estimate <- array(NA, dims)
  Estimate <- matrix(NA, nrow = nrow(X2), ncol = ncol(X2))
  # for(i in (1+p):dims[3]){
  for(i in (1+p):ncol(X2)){
    ## Eq 3
    # Estimate[,,i] <-
    #   (sigma_sqrt_list[[i]] %*% C_eps %*% sigma_sqrt_list[[i]]) /
    #   dims[1]^2
    Estimate[,i] <-
      as.vector( (sigma_sqrt_list[[i]] %*% C_eps %*% sigma_sqrt_list[[i]]) ) /
      nrow(X2)^2
  }

  Estimate
}

.predict_ARCH_op_or1 <- function(X2, alpha, delta, C_eps, p){
  # dims <- dim(X2)
  resolution <- sqrt(nrow(X2))

  # Eq 3 (Eq 29)
  # Estimate <- array(NA, dims)
  Estimate <- matrix(NA, nrow = nrow(X2), ncol = ncol(X2))
  sim_sigma2_matrix1 <- matrix(0,ncol = ncol(X2),nrow=resolution^2)
  # for(i in (1+p):dims[3]){
  for(i in (1+p):ncol(X2)){
    tmp_alpha1 <- rep(0, resolution^2)
    for(pp in 1:p){
      # prev_tensor <-
      #   as.vector( X2[,i-pp] %*% t(X2[,i-pp]) ) / resolution

      tmp_alpha1 <- tmp_alpha1 + alpha %*% X2[,i-pp] / resolution
    }

    ## Eq 3
    sim_sigma2_matrix1[,i] <- delta + tmp_alpha1

    tryCatch({
      Sigma2_svd1 <- svd(matrix(sim_sigma2_matrix1[,i],
                                ncol=resolution,nrow=resolution))
    }, error=function(e){
      cat(i)
      View(sim_sigma2_matrix1)
      View(arch_data1)
      stop('Error')
    })
    Sigma2_sqrt1 <- sqrt(resolution) * Sigma2_svd1$u %*%
      diag(sqrt(Sigma2_svd1$d)) %*% t(Sigma2_svd1$v)
    Estimate[,i] <- Sigma2_sqrt1 %*% C_eps %*% Sigma2_sqrt1 / resolution^2
  }

  Estimate
}
#########################################################################
#' Title
#'
#' @param resolution
#' @param N_fit
#' @param N_pred
#' @param alpha
#' @param delta
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
generateFitCompareOPARCH5 <-
  function(
    resolution = 20, N_fit = 500, N_pred = 250,
    alpha = function(s,t,u,v){
      val_seq <- seq(0,1,length.out=20)
      K <- 10 #30 #20 #15 #10 #5

      # tmp <- fda::create.fourier.basis(nbasis = K,rangeval = c(0,1))
      # basis_vals <- fda::eval.basis(seq(0,1,length.out=length(val_seq)), tmp)
      #
      # alp <- matrix(0,nrow=length(val_seq)^2,ncol = length(val_seq)^2)
      # for(k in 1:K){
      #   hkk <- basis_vals[,k] %*% t(basis_vals[,k]) /
      #     nrow(basis_vals)
      #   # tmp <- max(hkk)
      #   alp <- alp + kronecker(hkk, hkk) / nrow(hkk)
      # }
      # alp <- 2500*(alp / K)*20
      # saveRDS(alp,
      #         paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'_opt1.rds'))
      alp <- readRDS(
        paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'_opt1.rds'))
      alp[(which(s==val_seq)-1)*length(val_seq)+which(t==val_seq),
          (which(u==val_seq)-1)*length(val_seq)+which(v==val_seq)]
    },
    delta = function(s,t){pmin(t,s)+1 },
    seed=NA){
    if(!is.na(seed)) set.seed(seed)

    dat <- generateFunctionalCondVol_ff(resolution = resolution,
                                        N = N_fit+N_pred,'arch',
                                        alpha_par = alpha,
                                        delta_par = delta,
                                        p = 5)
    # dat_fit <- dat$garch_mat[,1:N_fit]
    # dat_pred <- dat$garch_mat[,N_fit+ ((1-5):N_pred)]
    dat_fit <- dat$data[,1:N_fit]
    dat_pred <- dat$data[,N_fit+ ((1-5):N_pred)]

    ######
    ## MLE
    ######
    NBASIS = 20
    fd = fda::Data2fd(argvals=seq(0,1,len=nrow(dat$data)),#len=nrow(dat$garch_mat)),
                      y=dat_fit,
                      fda::create.bspline.basis(nbasis=NBASIS))

    # extract data-driven basis functions through the truncated FPCA method.
    basis_est = basis.est(yd = dat_fit, M=2, type = "tfpca",NBASIS=NBASIS)$basis

    # estimate an FARCH(1) model with basis when M=1.
    arch1_est = est.fArch(fdata = fd, basis = basis_est[,1],q = 5)
    alpha_est1 <- as.numeric(arch1_est$As[[1]]) * (as.matrix(basis_est[,1]) %*% t(as.matrix(basis_est[,1])))
    alpha_est2 <- as.numeric(arch1_est$As[[2]]) * (as.matrix(basis_est[,1]) %*% t(as.matrix(basis_est[,1])))
    alpha_est3 <- as.numeric(arch1_est$As[[3]]) * (as.matrix(basis_est[,1]) %*% t(as.matrix(basis_est[,1])))
    alpha_est4 <- as.numeric(arch1_est$As[[4]]) * (as.matrix(basis_est[,1]) %*% t(as.matrix(basis_est[,1])))
    alpha_est5 <- as.numeric(arch1_est$As[[5]]) * (as.matrix(basis_est[,1]) %*% t(as.matrix(basis_est[,1])))
    delta_est <- basis_est[,1] * arch1_est$d

    result_mle <- list('hat_alpha'=
                         list(alpha_est1, alpha_est2,
                              alpha_est3, alpha_est4,
                              alpha_est5),
                       'hat_delta'=delta_est)

    ######
    ## Pointwise
    ######
    result_pw <- pointwiseARCHp(X = dat_fit, w = 0.01,p = 5)

    ######
    ## Operator - YW
    ######
    result_yw <- operatorARCH(X = dat_fit, w = 0.01,p = 5 )

    ######
    ## Operator - MP
    ######
    result_mp80 <- operatorARCHmp(X = dat_fit, w = 0.2, p = 5 )
    result_mp90 <- operatorARCHmp(X = dat_fit, w = 0.1, p = 5 )
    result_mp95 <- operatorARCHmp(X = dat_fit, w = 0.05,p = 5 )
    result_mp99 <- operatorARCHmp(X = dat_fit, w = 0.01,p = 5 )

    ######
    ## Get error for all
    ######

    # C_eps
    C_eps <- matrix(NA,resolution,resolution)
    times <- 1:resolution/resolution
    for (i in 1:resolution){
      #C_eps[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
      C_eps[i,] <- exp(-abs(times[i]-times)/2)
    }

    dat_pred2_pw <- dat_pred^2
    # dat_pred2 <- array(apply(dat_pred, 2, tcrossprod),
    #                    dim=c(nrow(dat_pred),nrow(dat_pred),ncol(dat_pred)))
    dat_pred2_old <- array(apply(dat_pred, 2, tcrossprod),
                        dim=c(nrow(dat_pred),
                              nrow(dat_pred),
                              ncol(dat_pred)))
    dat_pred2 <- matrix(nrow=nrow(dat_pred)^2,ncol=ncol(dat_pred)) # (nn x N)
    for(i in 1:ncol(dat_pred)){
      dat_pred2[,i] <- as.vector(dat_pred2_old[,,i])
    }
    estimate_or <-
      .predict_ARCH_op_or(X2 = dat_pred2,
                          sigma_sqrt_list=dat$params$Sigma2_sqrt_list[N_fit+(1-5):N_pred],
                          C_eps = C_eps, p = 5)
    estimate_parm <-
      .predict_ARCH_op_or1(X2 = dat_pred2,
                           alpha=dat$params$alpha,
                           delta=dat$params$delta,
                          C_eps = C_eps, p = 5)

    estimate_yw <-
      .predict_ARCH_op(X2=dat_pred2,
                       alpha=result_yw$hat_alpha,
                       delta=result_yw$hat_delta,
                       C_eps=C_eps, p=5)

    estimate_mp80 <-
      .predict_ARCH_op(X2=dat_pred2,
                       alpha=result_mp80$hat_alpha,
                       delta=result_mp80$hat_delta,
                       C_eps=C_eps, p=5)
    estimate_mp90 <-
      .predict_ARCH_op(X2=dat_pred2,
                       alpha=result_mp90$hat_alpha,
                       delta=result_mp90$hat_delta,
                       C_eps=C_eps, p=5)
    estimate_mp95 <-
      .predict_ARCH_op(X2=dat_pred2,
                       alpha=result_mp95$hat_alpha,
                       delta=result_mp95$hat_delta,
                       C_eps=C_eps, p=5)
    estimate_mp99 <-
      .predict_ARCH_op(X2=dat_pred2,
                       alpha=result_mp99$hat_alpha,
                       delta=result_mp99$hat_delta,
                       C_eps=C_eps, p=5)

    estimate_pw <-
      .pointwise_next_step_ests_op5(X2 = dat_pred2_pw,
                                    alpha = result_pw$hat_alpha,
                                    delta = result_pw$hat_delta,
                                    C_eps = C_eps, p = 5)
    estimate_mle <-
      .pointwise_next_step_ests_op5(X2 = dat_pred2_pw,
                                    alpha = result_mle$hat_alpha,
                                    delta = result_mle$hat_delta,
                                    C_eps = C_eps, p = 5)

    ## Get Errors
    # full_dat2 <- array(
    #   apply(dat$garch_mat, 2, tcrossprod),
    #   dim=c(nrow(dat$garch_mat),nrow(dat$garch_mat),ncol(dat$garch_mat)))
    full_dat2 <- apply(dat$data, 2, tcrossprod)
    # full_dat2 <- apply(dat$garch_mat, 2, tcrossprod)
    error <-
      data.frame('OR'=rep(0, ncol(dat_pred2)),
                 'Parm'=0,
                 'YW'=0,
                 'MP80'=0,'MP90'=0,'MP95'=0,'MP99'=0,
                 'PW'=0,
                 'MLE'=0,
                 'Mean'=0,
                 'Prev'=0,
                 'Zero'=0)
    for(i in (1+5):N_pred){
      # Errors for next step (learning upto last, no re-fit)
      # error[i,] <- c(
      #   dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_or[,,i])^2)),
      #   dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_yw[,,i])^2)),
      #   dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_mp80[,,i])^2)),
      #   dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_mp90[,,i])^2)),
      #   dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_mp95[,,i])^2)),
      #   dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_mp99[,,i])^2)),
      #   dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_pw[,,i])^2)),
      #   dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_mle[,,i])^2)),
      #   dot_integrate(dot_integrate_col((dat_pred2[,,i] -
      #                                      rowMeans(full_dat2[,,1:(N_fit+i-2)],dims=2))^2)),
      #   dot_integrate(dot_integrate_col((dat_pred2[,,i] - dat_pred2[,,i-1])^2)),
      #   dot_integrate(dot_integrate_col((dat_pred2[,,i] - 0)^2))
      # )
      error[i,] <- c(
        dot_integrate((dat_pred2[,i] - estimate_or[,i])^2),
        dot_integrate((dat_pred2[,i] - estimate_parm[,i])^2),
        dot_integrate((dat_pred2[,i] - estimate_yw[,i])^2),
        dot_integrate((dat_pred2[,i] - estimate_mp80[,i])^2),
        dot_integrate((dat_pred2[,i] - estimate_mp90[,i])^2),
        dot_integrate((dat_pred2[,i] - estimate_mp95[,i])^2),
        dot_integrate((dat_pred2[,i] - estimate_mp99[,i])^2),
        dot_integrate((dat_pred2[,i] - estimate_pw[,i])^2),
        dot_integrate((dat_pred2[,i] - estimate_mle[,i])^2),
        dot_integrate((dat_pred2[,i] - rowMeans(full_dat2[,1:(N_fit+i-2)]))^2),
        dot_integrate((dat_pred2[,i] - dat_pred2[,i-1])^2),
        dot_integrate((dat_pred2[,i] - 0)^2)
      )
    }

    colSums(error)
  }

.pointwise_next_step_ests_op5 <- function(X2, alpha, delta, C_eps,p){

  estimate <- .predictARCHp(X2, alpha, delta, p = p, est_sqr=T)


  # Estimate <- array(NA, dim=c(nrow(X2),nrow(X2),ncol(X2)))
  Estimate <- matrix(NA, nrow=nrow(X2)^2, ncol=ncol(X2) )
  for(i in 1:ncol(estimate) ){
    sqrt_est <- sqrt( pmax(estimate[,i], 0) )
    # Estimate[,,i] <- t(sqrt_est %*% C_eps) %*% sqrt_est / length(sqrt_est)
    Estimate[,i] <-
      as.vector( t(sqrt_est %*% C_eps) %*% sqrt_est ) /
      length(sqrt_est)
  }

  Estimate
}

# .predict_ARCH5_op_or <- function(X2, sigma_sqrt_list, C_eps, p){
#   dims <- dim(X2)
#
#   # Eq 3 (Eq 29)
#   Estimate <- array(NA, dims)
#   for(i in (1+p):dims[3]){
#
#     ## Eq 3
#     Estimate[,,i] <- array(0,dim=dims[1:2])
#     for(j in 1:p){
#       ## Eq 3
#       Estimate[,,i] <- Estimate[,,i] +
#         (sigma_sqrt_list[[i-j]] %*% C_eps %*% sigma_sqrt_list[[i-j]]) /
#         dims[1]^2
#     }
#   }
#
#   Estimate
# }
