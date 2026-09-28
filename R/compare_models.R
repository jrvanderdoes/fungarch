

compare_models <- function(N=750, resolution=20,
                           alpha_scale = c(0.7,0.7),
                           est_dimen=3,
                           # c_eps_type = c('Full','Middle','Simple'),
                           # gen_dimen = 3,
                           # simple_eigs = rep(1,gen_dimen),
                           N_fit = 500, w=0.1, silent=FALSE){
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
  # C_eps (use)
  C_eps_fit <- matrix(0, nrow=resolution, ncol = resolution)

  # if(c_eps_type=='Full'){
  #   C_eps_fit <- C_eps
  #   eigs_fit <- eigenvalues
  #   basis_gen <- basis
  #
  #   eigs_gen <- NULL
  # } else if(c_eps_type=='Middle'){
  #   for(k in 1:gen_dimen){
  #     C_eps_fit <- C_eps_fit +
  #       eigenvalues[k] * (basis[,k] %o% basis[,k])
  #   }
  #   eigs_fit <- eigen( C_eps_fit )$values / resolution
  #   basis_gen <- basis[,1:gen_dimen,drop=FALSE]
  #
  #   eigs_gen <- eigenvalues
  # } else if(c_eps_type=='Simple'){
  #   for(k in 1:gen_dimen){
  #     C_eps_fit <- C_eps_fit +
  #       gen_eigs[k] * (basis[,k] %o% basis[,k])
  #   }
  #   eigs_fit <- eigen( C_eps_fit )$values / resolution
  #   basis_gen <- basis[,1:gen_dimen,drop=FALSE]
  #
  #   eigs_gen <- gen_eigs
  # } else{
  #   stop('See c_eps_type.', call. = FALSE)
  # }

  ##################
  ## SETUP ALPHA/DELTA
  ##################
  #   Alpha is based on tensor sum of `alpha_dimen` basis functions
  alpha <- matrix(0, nrow=resolution^2, ncol = resolution^2)
  for(k in 1:length(alpha_scale)){
    hkk <- as.vector(basis[,k] %o% basis[,k])
    alpha <- alpha + alpha_scale[k] * (hkk %o% hkk)
  }

  #   Delta is full C_eps
  delta <- .generate_ceps(resolution)


  ##################
  ## Generate Data
  ##################
  data <- generateFARCH_use(N = N, resolution = resolution,
                            alpha=alpha, delta=delta,
                            basis=basis,
                            burnin=100, p=1,
                            eigs=eigenvalues[1:length(alpha_scale)],#eigs_gen,
                            C_eps=C_eps)

  ##################
  ## CV
  ##################
  X2 <- array(apply(data$data1, 2, tcrossprod),
              dim=c(nrow(data$data1),nrow(data$data1),ncol(data$data1)))
  errors <- errors1 <- matrix(NA, ncol=6,nrow=N-N_fit)
  colnames(errors) <- c('oracle','diag','pw', 'last','mean','zero')
  for(i in (N_fit+1):N){
    if(!silent) cat(i-N_fit,', ')
    X2_prev <- X2[,,i-1]
    ## Oracle model
    oracle_pred <- .pred_diag(resolution, alphas = alpha_scale, basis = basis,
                              delta = delta, X2_prev = X2_prev, C_eps = C_eps)

    ## Diag model with correct
    diag_mod <- .fit_diag_model(X = data$data1[,1:(i-1)],
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
    # diag_pred2 <- .pred_diag(resolution, alphas = diag_mod$alpha,
    #                          basis=basis, delta=diag_mod$delta2, X2_prev=X2_prev,
    #                          C_eps=C_eps)


    ## PW Model - MP
    pw_mod <- .fit_pw(X=data$data1[,1:(i-1)], w=w)
    pw_pred <- .pred_pw(X=data$data1[,1:(i-1)], pw_mod$alpha,
                        pw_mod$delta, C_eps)

    # ## PW Model - Tik
    # pw_mod1 <- pointwiseARCH(X=data$data1[,1:(i-1)],w = 0.001)
    # pw_pred1 <- .pred_pw(X=data$data1[,1:(i-1)], pw_mod1$hat_alpha,
    #                      pw_mod1$hat_delta, C_eps)

    ## Last Value
    last_pred <- X2[,,i-1]
    ## Mean
    mean_pred <- rowMeans(X2[,,1:(i-1)], dims = 2)
    ## Zero
    zero_pred <- matrix(0,resolution,resolution)


    # Compute Errors
    #   TODO:: Verify this
    errors[i-N_fit,] <- c(
      dot_integrate( (X2[,,i] - oracle_pred/N)^2 ),
      dot_integrate( (X2[,,i] - diag_pred/N)^2 ),
      # dot_integrate( (X2[,,i] - diag_pred1/N)^2 ),
      # dot_integrate( (X2[,,i] - diag_pred2/N)^2 ),
      dot_integrate( (X2[,,i] - pw_pred)^2 ),
      # dot_integrate( (X2[,,i] - pw_pred1)^2 ),
      dot_integrate( (X2[,,i] - last_pred)^2 ),
      dot_integrate( (X2[,,i] - mean_pred)^2 ),
      dot_integrate( (X2[,,i] - zero_pred)^2 )
    )
    errors1[i-N_fit,] <- c(
      dot_integrate( abs(X2[,,i] - oracle_pred/N) ),
      dot_integrate( abs(X2[,,i] - diag_pred/N) ),
      # dot_integrate( abs(X2[,,i] - diag_pred1/N) ),
      # dot_integrate( abs(X2[,,i] - diag_pred2/N) ),
      dot_integrate( abs(X2[,,i] - pw_pred) ),
      dot_integrate( abs(X2[,,i] - last_pred) ),
      dot_integrate( abs(X2[,,i] - mean_pred) ),
      dot_integrate( abs(X2[,,i] - zero_pred) )
    )
  }

  list(errors,errors1)
}
