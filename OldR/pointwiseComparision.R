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
#' generateFitComparePWARCH(N_fit = 500, N_pred = 250,
#'   seed=123456,
#'   alpha = function(t,s){ 15 * t * (1-t) * s * (1-s) },
#'   delta = function(t){ rep(0.5,length(t)) })
#'
#' generateFitComparePWARCH(N_fit = 500, N_pred = 250,
#'   alpha = function(t,s){ 12 * t * (1-t) * s * (1-s) },
#'   delta = function(t){ rep(0.5,length(t)) })
generateFitComparePWARCH <-
  function(
    resolution = 20, N_fit = 1000, N_pred = 250,
    alpha = function(t,s){ 15 * t * (1-t) * s * (1-s) },
    delta = function(t){ rep(0.01,length(t)) },
    seed=NA){
    if(!is.na(seed)) set.seed(seed)

    dat <- generateFunctionalCondVol(resolution = resolution,
                                     N = N_fit+N_pred,'arch',
                                     alpha_par = alpha,
                                     delta_par = delta)
    dat_fit <- dat$garch_mat[,1:N_fit]
    dat_pred <- dat$garch_mat[,N_fit+ (0:N_pred)]

    ######
    ## MLE
    ######
    fd = fda::Data2fd(argvals=seq(0,1,len=nrow(dat$garch_mat)),
                      y=dat_fit,
                      fda::create.bspline.basis(nbasis=32))

    # extract data-driven basis functions through the truncated FPCA method.
    basis_est = basis.est(dat_fit, M=2, "tfpca")$basis

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
    dat_pred2_op <- array(apply(dat_pred, 2, tcrossprod),
                       dim=c(nrow(dat_pred),nrow(dat_pred),ncol(dat_pred)))

    # tmp <- dat_pred - rowMeans(dat_pred)
    # X2_pred_op <- array(apply(tmp, 2, tcrossprod),
    #             dim=c(nrow(tmp),nrow(tmp),ncol(tmp)))

    ######
    ## Oracle
    ######
    result_or <- .get_oracle(dat)

    ######
    ## Get L2 error for all
    ######
    #dat_pred2 <- dat_pred^2

    estimate_or <-
      .pointwise_next_step_ests(dat_pred2_pw,
                                result_or$hat_alpha,
                                result_or$hat_delta)

    estimate_pw <-
      .pointwise_next_step_ests(dat_pred2_pw,
                                result_pw$hat_alpha,
                                result_pw$hat_delta)

    tmp <-
      .predict_ARCH_op(X2=dat_pred2_op,
                       alpha=result_yw$hat_alpha,
                       delta=result_yw$hat_delta,
                       C_eps=C_eps, p=1)
    estimate_yw <- matrix(NA,nrow(estimate_pw),ncol(estimate_pw))
    for(i in 1:ncol(estimate_yw)){
      estimate_yw[,i] <- diag(tmp[,,i])
    }

    tmp90 <-
      .predict_ARCH_op(X2=dat_pred2_op,
                       alpha=result_mp90$hat_alpha,
                       delta=result_mp90$hat_delta,
                       C_eps=C_eps, p=1)
    tmp95 <-
      .predict_ARCH_op(X2=dat_pred2_op,
                       alpha=result_mp95$hat_alpha,
                       delta=result_mp95$hat_delta,
                       C_eps=C_eps, p=1)
    tmp99 <-
      .predict_ARCH_op(X2=dat_pred2_op,
                       alpha=result_mp99$hat_alpha,
                       delta=result_mp99$hat_delta,
                       C_eps=C_eps, p=1)
    estimate_mp90 <- estimate_mp95 <- estimate_mp99 <-
      matrix(NA,nrow(estimate_pw),ncol(estimate_pw))
    for(i in 1:ncol(estimate_mp90)){
      estimate_mp90[,i] <- diag(tmp90[,,i])
      estimate_mp95[,i] <- diag(tmp95[,,i])
      estimate_mp99[,i] <- diag(tmp99[,,i])
    }


    # tmp <- .predict_ARCH_op(X2_pred_op,
    #                         alpha=result_op$hat_alpha,
    #                         delta=result_op$hat_delta,
    #                         C_eps=C_eps, p=1)
    # estimate_op <- matrix(NA,nrow(estimate_pw),ncol(estimate_pw))
    # for(i in 1:ncol(estimate_op)){
    #   estimate_op[,i] <- diag(tmp[,,i])
    # }

    estimate_mle <-
      .pointwise_next_step_ests(dat_pred2_pw,
                                result_mle$hat_alpha,
                                result_mle$hat_delta)
    # estimate_all_mle <-
    #   .pointwise_remaining_steps_ests(dat_pred2,
    #                                   result_mle$hat_alpha,
    #                                   result_mle$hat_delta)


    ## Get Errors
    error <- error_all <-
      data.frame('OR'=rep(0, ncol(dat_pred2_pw)),
                 'YW'=0,
                 'MP90'=0,
                 'MP95'=0,
                 'MP99'=0,
                 'PW'=0,
                 'MLE'=0,
                 'Mean'=0,
                 'Prev'=0,
                 'Zero'=0)
    for(i in 2:ncol(dat_pred2_pw)){
      # Errors for next step (learning upto last, no re-fit)
      error[i,] <- c(
        dot_integrate((dat_pred2_pw[,i] - estimate_or[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - estimate_yw[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - estimate_mp90[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - estimate_mp95[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - estimate_mp99[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - estimate_pw[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - estimate_mle[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - rowMeans(dat$garch_mat[,1:(N_fit+i-2)]^2))^2),
        dot_integrate((dat_pred2_pw[,i] - dat_pred2_pw[,i-1])^2),
        dot_integrate((dat_pred2_pw[,i] - 0)^2)
      )
    }
    colSums(error)
}


.get_oracle <- function(dat){
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

.pointwise_next_step_ests <- function(X2, alpha, delta){
  estimate <- matrix(NA, ncol = ncol(X2),nrow = nrow(X2))

  for(i in 2:ncol(X2)){
    # TODO:: Change to matrix mult
    for(j in 1:nrow(X2)){
      vector_op <- alpha[j, ] * X2[,(i-1)] # X^2
      estimate[j,i] <- delta[j] + dot_integrate(as.matrix(vector_op))
    }
  }

  estimate
}

.operator_next_step_ests <- function(X2, alpha, delta, C_eps){
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
    estimate[,i] <- diag(Estimate[,,i])
  }

  estimate
}

.pointwise_remaining_steps_ests <- function(X2, alpha, delta){
  estimate <- matrix(NA, ncol = ncol(X2),nrow = nrow(X2))
  estimate[,1] <- X2[,1]

  for(i in 2:ncol(X2)){
    # TODO:: Change to matrix mult
    for(j in 1:nrow(X2)){
      vector_op <- alpha[j, ] * estimate[,(i-1)] # X^2
      estimate[j,i] <- delta[j] + dot_integrate(as.matrix(vector_op))
    }
  }

  estimate
}

#######################################################
generateFitComparePWARCH5 <-
  function(
    resolution = 20, N_fit = 500, N_pred = 250,
    alpha = function(t,s){ 3 * t * (1-t) * s * (1-s) },
    delta = function(t){ rep(0.01,length(t)) },
    seed=NA){
    if(!is.na(seed)) set.seed(seed)

    dat <- generateFunctionalCondVol(resolution = resolution,
                                     N = N_fit+N_pred,'arch',
                                     alpha_par = alpha,
                                     delta_par = delta,
                                     p=5)
    dat_fit <- dat$garch_mat[,1:N_fit]
    dat_pred <- dat$garch_mat[,N_fit+ (0:N_pred)]

    ######
    ## MLE
    ######
    fd = fda::Data2fd(argvals=seq(0,1,len=nrow(dat$garch_mat)),
                      y=dat_fit,
                      fda::create.bspline.basis(nbasis=32))

    # extract data-driven basis functions through the truncated FPCA method.
    basis_est = basis.est(dat_fit, M=2, "tfpca")$basis

    # estimate an FARCH(1) model with basis when M=1.
    arch1_est = est.fArch(fd, basis_est[,1],q = 5)
    alpha_est1 <- as.numeric(arch1_est$As[[1]]) * (as.matrix(basis_est[,1]) %*% t(as.matrix(basis_est[,1])))
    alpha_est2 <- as.numeric(arch1_est$As[[2]]) * (as.matrix(basis_est[,1]) %*% t(as.matrix(basis_est[,1])))
    alpha_est3 <- as.numeric(arch1_est$As[[3]]) * (as.matrix(basis_est[,1]) %*% t(as.matrix(basis_est[,1])))
    alpha_est4 <- as.numeric(arch1_est$As[[4]]) * (as.matrix(basis_est[,1]) %*% t(as.matrix(basis_est[,1])))
    alpha_est5 <- as.numeric(arch1_est$As[[5]]) * (as.matrix(basis_est[,1]) %*% t(as.matrix(basis_est[,1])))
    delta_est <- basis_est[,1] * arch1_est$d

    result_mle <- list('hat_alpha1'=alpha_est1,
                       'hat_alpha2'=alpha_est2,
                       'hat_alpha3'=alpha_est3,
                       'hat_alpha4'=alpha_est4,
                       'hat_alpha5'=alpha_est5,
                       'hat_delta'=delta_est)

    ######
    ## Pointwise
    ######
    result_pw <- pointwiseARCHp(X = dat_fit, w = 0.01, p = 5)

    ######
    ## Operator - YW
    ######
    result_yw <- operatorARCH(X = dat_fit, w = 0.01,p = 5)

    ######
    ## Operator - MP
    ######
    result_mp90 <- operatorARCHmp(X = dat_fit, w = 0.1,p = 5)
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
    dat_pred2_op <- array(apply(dat_pred, 2, tcrossprod),
                          dim=c(nrow(dat_pred),nrow(dat_pred),ncol(dat_pred)))

    # tmp <- dat_pred - rowMeans(dat_pred)
    # X2_pred_op <- array(apply(tmp, 2, tcrossprod),
    #             dim=c(nrow(tmp),nrow(tmp),ncol(tmp)))

    ######
    ## Oracle
    ######
    result_or <- .get_oracle_archp(dat,p=5)

    ######
    ## Get L2 error for all
    ######
    #dat_pred2 <- dat_pred^2

    estimate_or <-
      .pointwise_next_step_ests_ARCHp(dat_pred2_pw,
                                      result_or$hat_alpha,
                                      result_or$hat_delta,
                                      p=5)

    estimate_pw <-
      .pointwise_next_step_ests_ARCHp(dat_pred2_pw,
                                      result_pw$hat_alpha,
                                      result_pw$hat_delta,
                                      p=5)


    tmp <-
      .predict_ARCH_op(X2=dat_pred2_op,
                       alpha=result_yw$hat_alpha,
                       delta=result_yw$hat_delta,
                       C_eps=C_eps, p=5)
    estimate_yw <- matrix(NA,nrow(estimate_pw),ncol(estimate_pw))
    for(i in 1:ncol(estimate_yw)){
      estimate_yw[,i] <- diag(tmp[,,i])
    }

    tmp90 <-
      .predict_ARCH_op(X2=dat_pred2_op,
                       alpha=result_mp90$hat_alpha,
                       delta=result_mp90$hat_delta,
                       C_eps=C_eps, p=5)
    tmp95 <-
      .predict_ARCH_op(X2=dat_pred2_op,
                       alpha=result_mp95$hat_alpha,
                       delta=result_mp95$hat_delta,
                       C_eps=C_eps, p=5)
    tmp99 <-
      .predict_ARCH_op(X2=dat_pred2_op,
                       alpha=result_mp99$hat_alpha,
                       delta=result_mp99$hat_delta,
                       C_eps=C_eps, p=5)
    estimate_mp90 <- estimate_mp95 <- estimate_mp99 <-
      matrix(NA,nrow(estimate_pw),ncol(estimate_pw))
    for(i in 1:ncol(estimate_mp90)){
      estimate_mp90[,i] <- diag(tmp90[,,i])
      estimate_mp95[,i] <- diag(tmp95[,,i])
      estimate_mp99[,i] <- diag(tmp99[,,i])
    }


    # tmp <- .predict_ARCH_op(X2_pred_op,
    #                         alpha=result_op$hat_alpha,
    #                         delta=result_op$hat_delta,
    #                         C_eps=C_eps, p=1)
    # estimate_op <- matrix(NA,nrow(estimate_pw),ncol(estimate_pw))
    # for(i in 1:ncol(estimate_op)){
    #   estimate_op[,i] <- diag(tmp[,,i])
    # }

    estimate_mle <-
      .pointwise_next_step_ests_ARCHp(dat_pred2_pw,
                                      result_mle$hat_alpha,
                                      result_mle$hat_delta,
                                      p=5)

    # estimate_all_mle <-
    #   .pointwise_remaining_steps_ests(dat_pred2,
    #                                   result_mle$hat_alpha,
    #                                   result_mle$hat_delta)


    ## Get Errors
    error <- error_all <-
      data.frame('OR'=rep(0, ncol(dat_pred2_pw)),
                 'YW'=0,
                 'MP90'=0,
                 'MP95'=0,
                 'MP99'=0,
                 'PW'=0,
                 'MLE'=0,
                 'Mean'=0,
                 'Prev'=0,
                 'Zero'=0)
    for(i in (5+1):ncol(dat_pred2_pw)){
      # Errors for next step (learning upto last, no re-fit)
      error[i,] <- c(
        dot_integrate((dat_pred2_pw[,i] - estimate_or[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - estimate_yw[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - estimate_mp90[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - estimate_mp95[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - estimate_mp99[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - estimate_pw[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - estimate_mle[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - rowMeans(dat$garch_mat[,1:(N_fit+i-2)]^2))^2),
        dot_integrate((dat_pred2_pw[,i] - dat_pred2_pw[,i-1])^2),
        dot_integrate((dat_pred2_pw[,i] - 0)^2)
      )
    }
    colSums(error)
  }


.get_oracle_archp <- function(dat,p){
  s <- seq(0,1,length.out=nrow(dat$garch_mat))
  alpha_true <- matrix(nrow=nrow(dat$garch_mat),ncol=nrow(dat$garch_mat))
  for(s1 in 1:length(s)){
    for(s2 in 1:length(s)){
      alpha_true[s1,s2] <- dat$params$alpha_par(s[s1],s[s2])
    }
  }

  d <- seq(0,1,length.out=nrow(dat$garch_mat))
  delta_true <- dat$params$delta_par(d)

  list('hat_alpha'=rep(list(alpha_true),p),
       'hat_delta'=delta_true)
}

.pointwise_next_step_ests_ARCHp <- function(X2, alpha, delta,p){
  estimate <- matrix(NA, ncol = ncol(X2),nrow = nrow(X2))

  for(i in (p+1):ncol(X2)){
    # TODO:: Change to matrix mult
    for(j in 1:nrow(X2)){
      estimate[j,i] <- delta[j]
      for(k in 1:p){
        vector_op <- alpha[[k]][j, ] * X2[,(i-p)] # X^2
        estimate[j,i] <- estimate[j,i] + dot_integrate(as.matrix(vector_op))
      }
    }
  }

  estimate
}

#######################################################

generateFitComparePWGARCH <-
  function(
    resolution = 20, N_fit = 1000, N_pred = 250,
    alpha = function(t,s){ 10 * t * (1-t) * s * (1-s) },
    beta = function(t,s){ 10 * t * (1-t) * s * (1-s) },
    delta = function(t){ rep(0.5,length(t)) },
    seed=NA){
    if(!is.na(seed)) set.seed(seed)

    dat <- generateFunctionalCondVol(resolution = resolution,
                                     N = N_fit+N_pred,'garch',
                                     alpha_par = alpha,
                                     beta_par = beta,
                                     delta_par = delta)
    dat_fit <- dat$garch_mat[,1:N_fit]
    dat_pred <- dat$garch_mat[,N_fit+ (0:N_pred)]

    ######
    ## MLE
    ######
    fd = fda::Data2fd(argvals=seq(0,1,len=nrow(dat$garch_mat)),
                      y=dat_fit,
                      fda::create.bspline.basis(nbasis=32))

    # extract data-driven basis functions through the truncated FPCA method.
    basis_est = basis.est(dat_fit, M=2, "tfpca")$basis

    # estimate an FARCH(1) model with basis when M=1.
    garch11_est = est.fGarch(fd, basis_est[,1])
    alpha_est <- as.numeric(garch11_est$As[[1]]) * (as.matrix(basis_est[,1]) %*% t(as.matrix(basis_est[,1])))
    beta_est <- as.numeric(garch11_est$Bs[[1]]) * (as.matrix(basis_est[,1]) %*% t(as.matrix(basis_est[,1])))
    delta_est <- basis_est[,1] * garch11_est$d

    result_mle <- list('hat_alpha'=alpha_est,
                       'hat_beta'=beta_est,
                       'hat_delta'=delta_est)

    ######
    ## Pointwise
    ######
    result_pw <- pointwiseGARCH(X = dat_fit, w = 0.01,L = 5)

    ######
    ## Operator - YW
    ######
    result_yw <- operatorGARCH(X = dat_fit, w = 0.01,L = 5)
    # result_yw1 <- operatorGARCH_2params(X = dat_fit, w = 0.01,L = 5)

    # ######
    # ## Operator - MP
    # ######
    # result_mp80 <- operatorGARCHmp(X = dat_fit, w = 0.2,L = 5)
    # result_mp90 <- operatorGARCHmp(X = dat_fit, w = 0.1,L = 5)
    # result_mp95 <- operatorGARCHmp(X = dat_fit, w = 0.05,L = 5)
    # result_mp99 <- operatorGARCHmp(X = dat_fit, w = 0.01,L = 5)

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
    dat_pred2_op <- array(apply(dat_pred, 2, tcrossprod),
                          dim=c(nrow(dat_pred),nrow(dat_pred),ncol(dat_pred)))

    # tmp <- dat_pred - rowMeans(dat_pred)
    # X2_pred_op <- array(apply(tmp, 2, tcrossprod),
    #             dim=c(nrow(tmp),nrow(tmp),ncol(tmp)))

    ######
    ## Oracle
    ######
    result_or <- .get_oracle_garch(dat)

    ######
    ## Get L2 error for all
    ######
    #dat_pred2 <- dat_pred^2

    estimate_or <-
      .pointwise_next_step_garch_ests(dat_pred2_pw,
                                      result_or$hat_alpha,
                                      result_or$hat_beta,
                                      result_or$hat_delta)

    estimate_pw <-
      .pointwise_next_step_garch_ests(dat_pred2_pw,
                                      result_pw$hat_alpha,
                                      result_pw$hat_beta,
                                      result_pw$hat_delta)

    tmp <-
      .predict_GARCH_op(X2=dat_pred2_op,
                        alpha=result_yw$hat_alpha,
                        beta=result_yw$hat_beta,
                        delta=result_yw$hat_delta,
                        C_eps=C_eps)
    estimate_yw <- matrix(NA,nrow(estimate_pw),ncol(estimate_pw))
    for(i in 1:ncol(estimate_yw)){
      estimate_yw[,i] <- diag(tmp[,,i])
    }

    # tmp80 <-
    #   .predict_GARCH_op(X2=dat_pred2_op,
    #                     alpha=result_mp80$hat_alpha,
    #                     beta=result_mp80$hat_beta,
    #                     delta=result_mp80$hat_delta,
    #                     C_eps=C_eps)
    #
    # tmp90 <-
    #   .predict_GARCH_op(X2=dat_pred2_op,
    #                     alpha=result_mp90$hat_alpha,
    #                     beta=result_mp90$hat_beta,
    #                     delta=result_mp90$hat_delta,
    #                     C_eps=C_eps)
    # tmp95 <-
    #   .predict_GARCH_op(X2=dat_pred2_op,
    #                     alpha=result_mp95$hat_alpha,
    #                     beta=result_mp95$hat_beta,
    #                     delta=result_mp95$hat_delta,
    #                     C_eps=C_eps)
    # tmp99 <-
    #   .predict_GARCH_op(X2=dat_pred2_op,
    #                     alpha=result_mp99$hat_alpha,
    #                     beta=result_mp99$hat_beta,
    #                     delta=result_mp99$hat_delta,
    #                     C_eps=C_eps)
    # estimate_mp80 <- estimate_mp90 <- estimate_mp95 <- estimate_mp99 <-
    #   matrix(NA,nrow(estimate_pw),ncol(estimate_pw))
    # for(i in 1:ncol(estimate_mp90)){
    #   estimate_mp80[,i] <- diag(tmp80[,,i])
    #   estimate_mp90[,i] <- diag(tmp90[,,i])
    #   estimate_mp95[,i] <- diag(tmp95[,,i])
    #   estimate_mp99[,i] <- diag(tmp99[,,i])
    # }


    # tmp <- .predict_ARCH_op(X2_pred_op,
    #                         alpha=result_op$hat_alpha,
    #                         delta=result_op$hat_delta,
    #                         C_eps=C_eps, p=1)
    # estimate_op <- matrix(NA,nrow(estimate_pw),ncol(estimate_pw))
    # for(i in 1:ncol(estimate_op)){
    #   estimate_op[,i] <- diag(tmp[,,i])
    # }

    estimate_mle <-
      .pointwise_next_step_garch_ests(dat_pred2_pw,
                                result_mle$hat_alpha,
                                result_mle$hat_beta,
                                result_mle$hat_delta)
    # estimate_all_mle <-
    #   .pointwise_remaining_steps_ests(dat_pred2,
    #                                   result_mle$hat_alpha,
    #                                   result_mle$hat_delta)


    ## Get Errors
    error <- error_all <-
      data.frame('OR'=rep(0, ncol(dat_pred2_pw)),
                 'YW'=0,
                 # 'MP80'=0,
                 # 'MP90'=0,
                 # 'MP95'=0,
                 # 'MP99'=0,
                 'PW'=0,
                 'MLE'=0,
                 'Mean'=0,
                 'Prev'=0,
                 'Zero'=0)
    for(i in 2:ncol(dat_pred2_pw)){
      # Errors for next step (learning upto last, no re-fit)
      error[i,] <- c(
        dot_integrate((dat_pred2_pw[,i] - estimate_or[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - estimate_yw[,i])^2),
        # dot_integrate((dat_pred2_pw[,i] - estimate_mp80[,i])^2),
        # dot_integrate((dat_pred2_pw[,i] - estimate_mp90[,i])^2),
        # dot_integrate((dat_pred2_pw[,i] - estimate_mp95[,i])^2),
        # dot_integrate((dat_pred2_pw[,i] - estimate_mp99[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - estimate_pw[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - estimate_mle[,i])^2),
        dot_integrate((dat_pred2_pw[,i] - rowMeans(dat$garch_mat[,1:(N_fit+i-2)]^2))^2),
        dot_integrate((dat_pred2_pw[,i] - dat_pred2_pw[,i-1])^2),
        dot_integrate((dat_pred2_pw[,i] - 0)^2)
      )
    }
    colSums(error)
  }


.get_oracle_garch <- function(dat){
  s <- seq(0,1,length.out=nrow(dat$garch_mat))
  alpha_true <- beta_true <-
    matrix(nrow=nrow(dat$garch_mat),ncol=nrow(dat$garch_mat))
  for(s1 in 1:length(s)){
    for(s2 in 1:length(s)){
      alpha_true[s1,s2] <- dat$params$alpha_par(s[s1],s[s2])
    }
  }
  for(s1 in 1:length(s)){
    for(s2 in 1:length(s)){
      beta_true[s1,s2] <- dat$params$beta_par(s[s1],s[s2])
    }
  }

  d <- seq(0,1,length.out=nrow(dat$garch_mat))
  delta_true <- dat$params$delta_par(d)

  list('hat_delta'=delta_true,
       'hat_alpha'=alpha_true,
       'hat_beta'=beta_true)
}

.pointwise_next_step_garch_ests <- function(X2, alpha, beta, delta){
  sigma2 <- matrix(NA, ncol = ncol(X2),nrow = nrow(X2))

  sigma2[,1] <- 0
  for(i in 2:ncol(X2)){
    # TODO:: Change to matrix mult
    for(j in 1:nrow(X2)){
      vector_op <- alpha[j, ] * X2[,(i-1)] # X^2
      # vector_for_beta_op <- beta[j,] * X2[,(i-1)]
      vector_for_beta_op <- beta[j,] * sigma2[,(i-1)]
      sigma2[j,i] <- delta[j] +
        dot_integrate(as.matrix(vector_op)) +
        dot_integrate(as.matrix(vector_for_beta_op))
    }
  }

  sigma2
}

.operator_next_step_garch_ests <- function(X2, alpha, beta, delta, C_eps){
  dims <- dim(X2)
  Sigma <- array(NA, dims)
  Sigma[,,1] <- delta
  estimate <- data.frame(matrix(nrow=dims[2],ncol = dims[3]))

  for(i in 2:dims[3]){
    for(j in 1:dims[1]){
      vector_op <- alpha[j, ] * X2[,,(i-1)] # X^2
      vector_op_garch <- beta[j, ] * X2[,,(i-1)] # X^2
      Sigma[j,,i] <- delta[j,] +
        colSums(vector_op)/dims[1] +
        colSums(vector_op_garch)/dims[1]
    }
  }

  Estimate <- array(NA, dims)
  for(i in 1:dims[3]){
    Sigma_svd <- svd(Sigma[,,i])
    Sigma_sqrt <- Sigma_svd$u %*% diag(sqrt(Sigma_svd$d)) %*% t(Sigma_svd$v)

    Estimate[,,i] <- Sigma_sqrt %*% (C_eps) %*% Sigma_sqrt
    estimate[,i] <- diag(Estimate[,,i])
  }

  estimate
}
