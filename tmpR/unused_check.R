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
generateFitCompareOPARCH1 <-
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
    ## Operator
    ######
    result_op <- operatorARCH(X = dat_fit, w = 0.01,p = 2 )


    ######
    ## Get error for all
    ######

    # C_eps
    C_eps <- matrix(NA,resolution,resolution)
    times <- 1:resolution/resolution
    for (i in 1:resolution){
      C_eps[i,] <- exp(-times[i]/2-times/2) *
        pmin(exp(times[i]), exp(times))
    }

    dat_pred2_pw <- dat_pred^2
    dat_pred2 <- array(apply(dat_pred, 2, tcrossprod),
                       dim=c(nrow(dat_pred),nrow(dat_pred),ncol(dat_pred)))
    estimate_or <-
      .predict_ARCH_op_or(X2 = dat_pred2,
                          sigma_sqrt_list=dat$params$Sigma2_sqrt_list[N_fit+0:N_pred],
                          C_eps = C_eps, p = 1)

    estimate_op <-
      .predict_ARCH_op(X2=dat_pred2,
                       alpha=result_op$hat_alpha,
                       delta=result_op$hat_delta,
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
                 'OP'=0,
                 'PW'=0,
                 'MLE'=0,
                 'Mean'=0,
                 'Prev'=0,
                 'Zero'=0)
    for(i in 2:dim(dat_pred2)[3]){

      # Errors for next step (learning upto last, no re-fit)
      error[i,] <- c(
        dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_or[,,i])^2)),
        dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_mle[,,i])^2)),
        dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_op[,,i])^2)),
        dot_integrate(dot_integrate_col((dat_pred2[,,i] - estimate_pw[,,i])^2)),
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
  dims <- dim(X2)

  # Eq 3 (Eq 29)
  Estimate <- array(NA, dims)
  for(i in (1+p):dims[3]){

    ## Eq 3
    Estimate[,,i] <-
      (sigma_sqrt_list[[i]] %*% C_eps %*% sigma_sqrt_list[[i]]) /
      dims[1]^2
  }

  Estimate
}

