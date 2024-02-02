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
    resolution = 50, N_fit = 1000, N_pred = 250,
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
    ## Operator
    ######
    result_op <- operatorARCH(X = dat_fit, w = 0.01)

    # C_eps
    C_eps <- matrix(NA,resolution,resolution)
    times <- 1:resolution/resolution
    for (i in 1:resolution){
      C_eps[i,] <- exp(-times[i]/2-times/2) *
        pmin(exp(times[i]), exp(times))
    }

    tmp <- dat_pred - rowMeans(dat_pred)
    X2_pred_op <- array(apply(tmp, 2, tcrossprod),
                dim=c(nrow(tmp),nrow(tmp),ncol(tmp)))

    # .computeARCHError_tensor_new(
    #   X = X2_pred_op,
    #   alpha = result_op$hat_alpha,
    #   delta = result_op$hat_delta,
    #   C_eps = C_eps)
    ######
    ## Oracle
    ######
    result_or <- .get_oracle(dat)

    ######
    ## Get L2 error for all
    ######
    dat_pred2 <- dat_pred^2
    estimate_or <-
      .pointwise_next_step_ests(dat_pred2,
                                result_or$hat_alpha,
                                result_or$hat_delta)
    estimate_pw <-
      .pointwise_next_step_ests(dat_pred2,
                                result_pw$hat_alpha,
                                result_pw$hat_delta)
    estimate_op <-
      .operator_next_step_ests(X2_pred_op,
                               result_op$hat_alpha,
                               result_op$hat_delta,
                               C_eps)
    estimate_mle <-
      .pointwise_next_step_ests(dat_pred2,
                                result_mle$hat_alpha,
                                result_mle$hat_delta)
    # estimate_all_mle <-
    #   .pointwise_remaining_steps_ests(dat_pred2,
    #                                   result_mle$hat_alpha,
    #                                   result_mle$hat_delta)


    ## Get Errors
    error <- error_all <-
      data.frame('OR'=rep(0, ncol(dat_pred2)),
                 'PW'=0,
                 'OP'=0,
                 'MLE'=0,
                 'Mean'=0,
                 'Prev'=0)
    for(i in 2:ncol(dat_pred2)){
      # Errors for next step (learning upto last, no re-fit)
      error[i,] <- c(
        dot_integrate(abs(dat_pred2[,i] - estimate_or[,i])),
        dot_integrate(abs(dat_pred2[,i] - estimate_pw[,i])),
        dot_integrate(abs(dat_pred2[,i] - estimate_op[,i])),
        dot_integrate(abs(dat_pred2[,i] - estimate_mle[,i])),
        dot_integrate(abs(dat_pred2[,i] - rowMeans(dat$garch_mat[,1:(N_fit+i-2)]^2))),
        dot_integrate(abs(dat_pred2[,i] - dat_pred2[,i-1]))
      )
      # # Errors for all pred learning only from fit
      # error_all[i,] <- c(
      #   dot_integrate(abs(dat_pred2[,i] - estimate_all_or[,i])),
      #   dot_integrate(abs(dat_pred2[,i] - estimate_all_pw[,i])),
      #   dot_integrate(abs(dat_pred2[,i] - estimate_all_mle[,i])),
      #   dot_integrate(abs(dat_pred2[,i] - rowMeans(dat_fit^2))),
      #   dot_integrate(abs(dat_pred2[,i] - dat_pred2[,1]))
      # )
    }

    # results <- rbind(colSums(error),
    #                  colSums(error_all))
    # rownames(results) <- c('OneStep','FromLastFit')

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
