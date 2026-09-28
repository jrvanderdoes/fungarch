compare_models_fit_each <- function(data0, N_fit=750,
                           est_dimen=3,
                           w=0.1, silent=FALSE){
  resolution <- nrow(data0)
  N <- ncol(data0)
  ##################
  ##  C_eps
  ##################
  ## Full C_eps
  C_eps <- .generate_ceps(resolution)
  eig_Ce <- eigen( C_eps )
  basis <- eig_Ce$vectors * sqrt(resolution)
  basis <- fda::eval.basis(seq(0,1,length.out=resolution),
                           fda::create.fourier.basis(nbasis=resolution) )

  eigenvalues <- eig_Ce$values / resolution

  ##################
  ## CV
  ##################
  X2 <- array(apply(data0, 2, tcrossprod),
              dim=c(nrow(data0),nrow(data0),ncol(data0)))
  errors <- errors1 <- matrix(NA, ncol=5,nrow=N-N_fit)
  colnames(errors) <- c('diag','pw', 'last','mean','zero')
  for(i in (N_fit+1):N){
    if(!silent) cat(i-N_fit,', ')
    X2_prev <- X2[,,i-1]

    ## Diag model with correct
    diag_mod <- .fit_diag_model(X = data0[,1:(i-1)],
                                C_eps = C_eps,
                                basis = basis,
                                eigenvalues = eigenvalues,
                                w=w, p=1, dimen=est_dimen)
    # TODO:: Delta Information
    diag_pred <- .pred_diag(resolution, alphas = diag_mod$alpha,
                            basis=basis, delta=diag_mod$delta, X2_prev=X2_prev,
                            C_eps=C_eps)


    ## PW Model - MP
    pw_mod <- .fit_pw(X=data0[,1:(i-1)], w=w)
    pw_pred <- .pred_pw(X=data0[,1:(i-1)], pw_mod$alpha,
                        pw_mod$delta, C_eps)

    ## Last Value
    last_pred <- X2[,,i-1]
    ## Mean
    mean_pred <- rowMeans(X2[,,1:(i-1)], dims = 2)
    ## Zero
    zero_pred <- matrix(0,resolution,resolution)


    # Compute Errors
    #   TODO:: Verify this
    errors[i-N_fit,] <- c(
      dot_integrate( (X2[,,i] - diag_pred/N)^2 ),
      dot_integrate( (X2[,,i] - pw_pred)^2 ),
      dot_integrate( (X2[,,i] - last_pred)^2 ),
      dot_integrate( (X2[,,i] - mean_pred)^2 ),
      dot_integrate( (X2[,,i] - zero_pred)^2 )
    )
    errors1[i-N_fit,] <- c(
      dot_integrate( abs(X2[,,i] - diag_pred/N) ),
      dot_integrate( abs(X2[,,i] - pw_pred) ),
      dot_integrate( abs(X2[,,i] - last_pred) ),
      dot_integrate( abs(X2[,,i] - mean_pred) ),
      dot_integrate( abs(X2[,,i] - zero_pred) )
    )
  }

  list(errors,errors1)
}

compare_models_fit_once <- function(data0, N_fit=750,
                                    est_dimen=3,
                                    w=0.1, silent=FALSE){
  resolution <- nrow(data0)
  N <- ncol(data0)
  ##################
  ##  C_eps
  ##################
  ## Full C_eps
  C_eps <- .generate_ceps(resolution)
  eig_Ce <- eigen( C_eps )
  basis <- eig_Ce$vectors * sqrt(resolution)
  basis <- fda::eval.basis(seq(0,1,length.out=resolution),
                           fda::create.fourier.basis(nbasis=resolution) )

  eigenvalues <- eig_Ce$values / resolution


  ##################
  ## data0
  ##################
  X2 <- array(apply(data0, 2, tcrossprod),
              dim=c(nrow(data0),nrow(data0),ncol(data0)))

  ##################
  ## Fit Model
  ##################
  ## Diag model with correct
  diag_mod <- .fit_diag_model(X = data0[,1:N_fit],
                              C_eps = C_eps,
                              basis = basis,
                              eigenvalues = eigenvalues,
                              w=w, p=1, dimen=est_dimen)
  pw_mod <- .fit_pw(X=data0[,1:N_fit], w=w)

  ## Last Value
  last_pred0 <- X2[,,N_fit]
  ## Mean
  mean_pred0<- rowMeans(X2[,,1:N_fit], dims = 2)
  ## Zero
  zero_pred0 <- matrix(0,resolution,resolution)

  ##################
  ## CV
  ##################
  errors <- errors1 <- matrix(NA, ncol=8,nrow=N-N_fit)
  colnames(errors) <- c('diag','pw', 'last0','last','mean0','mean','zero0','zero')
  for(i in (N_fit+1):N){
    if(!silent) cat(i-N_fit,', ')
    X2_prev <- X2[,,i-1]

    # TODO:: Delta Information
    diag_pred <- .pred_diag(resolution, alphas = diag_mod$alpha,
                            basis=basis, delta=diag_mod$delta, X2_prev=X2_prev,
                            C_eps=C_eps)


    ## PW Model - MP
    pw_pred <- .pred_pw(X=data0[,1:(i-1)], pw_mod$alpha,
                        pw_mod$delta, C_eps)

    ## Last Value
    last_pred <- X2[,,i-1]
    ## Mean
    mean_pred <- rowMeans(X2[,,1:(i-1)], dims = 2)
    ## Zero
    zero_pred <- matrix(0,resolution,resolution)


    # Compute Errors
    #   TODO:: Verify this
    errors[i-N_fit,] <- c(
      dot_integrate( (X2[,,i] - diag_pred/N)^2 ),
      dot_integrate( (X2[,,i] - pw_pred)^2 ),
      dot_integrate( (X2[,,i] - last_pred0)^2 ),
      dot_integrate( (X2[,,i] - last_pred)^2 ),
      dot_integrate( (X2[,,i] - mean_pred0)^2 ),
      dot_integrate( (X2[,,i] - mean_pred)^2 ),
      dot_integrate( (X2[,,i] - zero_pred0)^2 ),
      dot_integrate( (X2[,,i] - zero_pred)^2 )
    )
    errors1[i-N_fit,] <- c(
      dot_integrate( abs(X2[,,i] - diag_pred/N) ),
      dot_integrate( abs(X2[,,i] - pw_pred) ),
      dot_integrate( abs(X2[,,i] - last_pred0) ),
      dot_integrate( abs(X2[,,i] - last_pred) ),
      dot_integrate( abs(X2[,,i] - mean_pred0) ),
      dot_integrate( abs(X2[,,i] - mean_pred) ),
      dot_integrate( abs(X2[,,i] - zero_pred0) ),
      dot_integrate( abs(X2[,,i] - zero_pred) )
    )
  }

  list(errors,errors1)
}
