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
    resolution = 25, N_fit = 500, N_pred = 250,
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
    arch1_est = est.fArch(fd, basis_est[,1])
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
    result_op <- operatorARCH(X = dat_fit, w = 0.01)

    # tmp <- dat_pred - rowMeans(dat_pred)
    # X2_pred_op <- array(apply(tmp, 2, tcrossprod),
    #                     dim=c(nrow(tmp),nrow(tmp),ncol(tmp)))

    # .computeARCHError_tensor_new(
    #   X = X2_pred_op,
    #   alpha = result_op$hat_alpha,
    #   delta = result_op$hat_delta,
    #   C_eps = C_eps)
    ######
    ## Oracle
    ######

    ## TODO:: Change this
    #result_or <- .get_oracle_op(dat)

    ######
    ## Get L2 error for all
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
      .predict_next_step_oracle_op(dat_pred2,
                                dat$params$alpha_par,
                                dat$params$delta_par)
    estimate_pw <-
      .pointwise_next_step_ests_op(dat_pred2_pw,
                                result_pw$hat_alpha,
                                result_pw$hat_delta,
                                C_eps)
    estimate_op <-
      .operator_next_step_ests_op(dat_pred2,
                               result_op$hat_alpha,
                               result_op$hat_delta,
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
                 'PW'=0,
                 'OP'=0,
                 'MLE'=0,
                 'Mean'=0,
                 'Prev'=0)
    for(i in 2:dim(dat_pred2)[3]){

      # Errors for next step (learning upto last, no re-fit)
      error[i,] <- c(
        dot_integrate(abs(dat_pred2[,,i] - estimate_or[,i])),
        dot_integrate(abs(dat_pred2[,,i] - estimate_pw[,,i])),
        dot_integrate(abs(dat_pred2[,,i] - estimate_op[,,i])),
        dot_integrate(abs(dat_pred2[,,i] - estimate_mle[,,i])),
        dot_integrate(abs(dat_pred2[,,i] - rowMeans(full_dat2[,,1:(N_fit+i-2)],dims=2))),
        dot_integrate(abs(dat_pred2[,,i] - dat_pred2[,,i-1]))
      )
    }

    # results <- rbind(colSums(error),
    #                  colSums(error_all))
    # rownames(results) <- c('OneStep','FromLastFit')

    colSums(error)
}


.get_oracle_op <- function(dat){
  s <- seq(0,1,length.out=nrow(dat$garch_mat))
  alpha_true <- matrix(nrow=nrow(dat$garch_mat),ncol=nrow(dat$garch_mat))
  for(s1 in 1:length(s)){
    for(s2 in 1:length(s)){
      alpha_true[s1,s2] <- dat$params$alpha_par(s[s1],s[s2])
    }
  }

  d <- seq(0,1,length.out=nrow(dat$garch_mat))
  delta_true <- dat$params$delta_par(d)

  list('hat_delta'=delta_true,
       'hat_alpha'=alpha_true)
}

.pointwise_next_step_ests_op <- function(X2, alpha, delta, C_eps){
  dims <- dim(X2)
  estimate <- data.frame(matrix(ncol=ncol(X2),nrow=nrow(X2)))
  estim <- array(NA, dim=c(nrow(X2),nrow(X2),ncol(X2)))

  for(i in 2:ncol(X2)){
    # TODO:: Change to matrix mult
    for(j in 1:nrow(X2)){
      vector_op <- alpha[j, ] * X2[,(i-1)] # X^2
      estimate[j,i] <- delta[j] + dot_integrate(as.matrix(vector_op))
    }
    estim[,,i] <- ( estimate[,i] %*% t(estimate[,i]) ) * C_eps
  }

  estim
}

.operator_next_step_ests_op <- function(X2, alpha, delta, C_eps){
  dims <- dim(X2)
  Sigma <- array(NA, dims)
  Sigma[,,1] <- delta
  estimate <- data.frame(matrix(nrow=dims[2],ncol = dims[3]))

  for(i in 2:dims[3]){
    for(j in 1:dims[1]){
      vector_op <- alpha[j, ] * X2[,,(i-1)] # X^2
      Sigma[j,,i] <- delta[j,] + colSums(vector_op)/dims[1]
    }
  }

  Estimate <- array(NA, dims)
  for(i in 1:dims[3]){
    Sigma_svd <- svd(Sigma[,,i])
    Sigma_sqrt <- Sigma_svd$u %*% diag(sqrt(Sigma_svd$d)) %*% t(Sigma_svd$v)

    Estimate[,,i] <- Sigma_sqrt %*% (C_eps) %*% Sigma_sqrt
  }

  Estimate
}
