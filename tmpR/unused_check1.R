
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
    result_op <- operatorARCH(X = dat_fit, w = 0.01,p = 1 )


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
        dot_integrate(dot_integrate_col(abs(dat_pred2[,,i] - estimate_or[,,i]))),
        dot_integrate(dot_integrate_col(abs(dat_pred2[,,i] - estimate_mle[,,i]))),
        dot_integrate(dot_integrate_col(abs(dat_pred2[,,i] - estimate_op[,,i]))),
        dot_integrate(dot_integrate_col(abs(dat_pred2[,,i] - estimate_pw[,,i]))),
        dot_integrate(dot_integrate_col(abs(dat_pred2[,,i] -
                                           rowMeans(full_dat2[,,1:(N_fit+i-2)],dims=2)))),
        dot_integrate(dot_integrate_col(abs(dat_pred2[,,i] - dat_pred2[,,i-1]))),
        dot_integrate(dot_integrate_col(abs(dat_pred2[,,i] - 0)))
      )
    }

    colSums(error)
  }



nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'OP'=NA,
             'PW'=NA,
             'MLE'=NA,
             'Mean'=NA,'Prev'=NA,'Zero'=NA)


cat('Fit (',nSims,'): ')
set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitCompareOPARCH1(N_fit = 500, N_pred = 250,
                             alpha = function(s,t,u,v){
                               f1 <- function(x){
                                 # int (x(1-x))^2 dx = 1/30
                                 x *(1-x) * 30
                               }
                               1/400 *
                                 (f1(s)*f1(t)*f1(u)*f1(v) ) },
                             delta = function(s,t){10*s*(1-s)*t*(1-t) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}

round(colMeans(results_oneStep,na.rm = T),2)
apply(results_oneStep,MARGIN = 2,
      function(X){round(median(X),3)})
