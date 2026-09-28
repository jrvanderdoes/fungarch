generateFARCH_use <- function(N, resolution, alpha, delta, basis,
                          burnin=100, p=1, eigs=NULL,
                          C_eps=NULL){
  # #### Functions
  # .get_kl_errors <- function(N, resolution, basis, eigs){
  #
  #   .get_errors <- function(basis,eigs){
  #     ## NEXT
  #     xi <- sapply(1:ncol(basis),
  #                  function(e,eigs) {
  #                    stats::rnorm(1, mean = 0, sd = sqrt(eigs[e]))
  #                  }, eigs=eigs
  #     )
  #
  #     Zeta <- tryCatch(basis * xi[col(basis)],
  #                      error = function(e) {
  #                        stop(call. = F, paste0(
  #                          "Check number of eigenvalues given. ",
  #                          "It does not match number of basis functions. ",
  #                          "Note, did you account for the constant function if ",
  #                          "it is in the basis?"
  #                        ))
  #                      }
  #     )
  #
  #     rowSums(Zeta)
  #   }
  #
  #   data <- data.frame(matrix(NA, ncol = N, nrow = resolution))
  #   for (j in 1:N) {
  #     data[, j] <- .get_errors(basis,eigs)
  #   }
  #
  #   data
  # }
  #### Functions
  .get_kl_errors <- function(N, resolution, basis, eigs){

    .get_errors <- function(basis,eigs){
      ## NEXT
      xi <- sapply(eigs,
                   function(e) {
                     stats::rnorm(1, mean = 0, sd = sqrt(e))
                   }
      )

      Zeta <- tryCatch(basis[,1:length(eigs),drop=FALSE] *
                         xi[col(basis[,1:length(eigs),drop=FALSE])],
                       error = function(e) {
                         stop(call. = F, paste0(
                           "Check number of eigenvalues given. ",
                           "It does not match number of basis functions. ",
                           "Note, did you account for the constant function if ",
                           "it is in the basis?"
                         ))
                       }
      )

      rowSums(Zeta)
    }

    data <- data.frame(matrix(NA, ncol = N, nrow = resolution))
    for (j in 1:N) {
      data[, j] <- .get_errors(basis,eigs)
    }

    data
  }




  #### Code
  times <- seq(0, 1, length.out=resolution)

  NBurn <- burnin + N
  sim_sigma2_matrix1 <-
    sim_sigma2_matrix2 <-
    matrix(0, ncol = NBurn, nrow = resolution^2)
  arch_data1 <- arch_data2 <-
    matrix(NA, nrow = resolution, ncol = NBurn)

  # if(is.null(eigs)){
    # eps <- t(MASS::mvrnorm(n = NBurn,
    #                        mu = c(rep(0,resolution)),
    #                        Sigma = C_eps,
    #                        empirical = FALSE))
  # }else{
    eps <- .get_kl_errors(N = NBurn, resolution = resolution,
                          basis=basis, eigs=eigs)
  # }

  ## Define variance and GARCH matrices
  # Initialize
  for(i in 1:p){
    # arch_data1[,i] <- operator(delta, eps[,i])
    arch_data2[,i] <- operator(delta, eps[,i])
  }


  # Set matrices
  Sigma2_list1 <- Sigma2_sqrt_list1 <-
    Sigma2_list2 <- Sigma2_sqrt_list2 <- list()

  # Compute values
  for(i in (p+1):NBurn){
    tmp_alpha1 <- tmp_alpha2 <- rep(0, resolution^2)
    for(pp in 1:p){
      # # Alpha 1
      # prev_tensor <- tensor_vec(arch_data1[,i-pp], arch_data1[,i-pp])
      # for(b in 2:K){
      #   mat_val <- basis[,b] %o% basis[,b]
      #   tmp <- mat_val * scale *
      #     as.numeric((basis[,b] %*% arch_data1[,i-pp] / nrow(basis))^2)
      #   tmp_alpha1 <- tmp_alpha1 + tmp
      # }

      # Alpha 2
      prev_tensor <- tensor_vec(arch_data2[,i-pp], arch_data2[,i-pp])
      tmp_alpha2 <- tmp_alpha2 + operator(alpha, prev_tensor)
    }

    #  Sigma^2
    # sim_sigma2_matrix1[,i] <- as.vector(delta) + tmp_alpha1
    sim_sigma2_matrix2[,i] <- as.vector(delta) + tmp_alpha2

    # Compute the sqrt matrix
    tryCatch({
      # Sigma2_svd1 <- svd(matrix(sim_sigma2_matrix1[,i],
      #                           ncol=resolution,nrow=resolution))
      Sigma2_svd2 <- svd(matrix(sim_sigma2_matrix2[,i],
                                ncol=resolution,nrow=resolution))
    }, error=function(e){
      stop('Error in sqrt matrix')
    })
    Sigma2_sqrt2 <- sqrt(resolution) *
      Sigma2_svd2$u %*% diag(sqrt(Sigma2_svd2$d)) %*% t(Sigma2_svd2$v)
    arch_data2[,i] <- operator(Sigma2_sqrt2, eps[,i])

    # Sigma2_sqrt_list1[[i]] <- Sigma2_sqrt1
    Sigma2_sqrt_list2[[i]] <- Sigma2_sqrt2
  }

  # Return
  new_res1 <- list(#data = arch_data1[, burnin + 1:N],
    data1 = arch_data2[, burnin + 1:N],
    #sigma_mat = sim_sigma2_matrix1[, burnin + 1:N],
    sigma_mat1 = sim_sigma2_matrix2[, burnin + 1:N],
    error_matrix = eps[, burnin + 1:N],
    params=list(p=p, alpha = alpha, delta = delta)
  )
}

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
  errors <- matrix(NA, ncol=9,nrow=N-N_fit)
  colnames(errors) <- c('oracle','diag','diag1','diag2','pw', 'pw-t','last','mean','zero')
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
    diag_pred1 <- .pred_diag(resolution, alphas = diag_mod$alpha,
                            basis=basis, delta=diag_mod$delta1, X2_prev=X2_prev,
                            C_eps=C_eps)
    diag_pred2 <- .pred_diag(resolution, alphas = diag_mod$alpha,
                            basis=basis, delta=diag_mod$delta2, X2_prev=X2_prev,
                            C_eps=C_eps)


    ## PW Model - MP
    pw_mod <- .fit_pw(X=data$data1[,1:(i-1)], w=w)
    pw_pred <- .pred_pw(X=data$data1[,1:(i-1)], pw_mod$alpha,
                        pw_mod$delta, C_eps)

    ## PW Model - Tik
    pw_mod1 <- pointwiseARCH(X=data$data1[,1:(i-1)],w = 0.001)
    pw_pred1 <- .pred_pw(X=data$data1[,1:(i-1)], pw_mod1$hat_alpha,
                         pw_mod1$hat_delta, C_eps)

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
      dot_integrate( (X2[,,i] - diag_pred1/N)^2 ),
      dot_integrate( (X2[,,i] - diag_pred2/N)^2 ),
      dot_integrate( (X2[,,i] - pw_pred)^2 ),
      dot_integrate( (X2[,,i] - pw_pred1)^2 ),
      dot_integrate( (X2[,,i] - last_pred)^2 ),
      dot_integrate( (X2[,,i] - mean_pred)^2 ),
      dot_integrate( (X2[,,i] - zero_pred)^2 )
    )
  }

  errors
}


###################################
# res <- compare_models(N=550)
# colMeans(res)
C_eps <- .generate_ceps(20)
eigenvalues <- eigen( C_eps )$values / 20

nSims <- 50
results <- matrix(NA, ncol=9,nrow=nSims)
colnames(results) <- c('oracle','diag','diag1','diag2','pw','pw1','last','mean','zero')
for(i in 1:nSims){
  cat(i, ', ')
  set.seed(1234+i)
  res <- tryCatch({
    compare_models(N = 1250,N_fit = 1000,silent=TRUE,
                   alpha_scale = c(rep(0.9,3))/eigenvalues[1:3])
  }, error=function(e){
    compare_models(N = 1250,N_fit = 1000,silent=TRUE,
                   alpha_scale =rep(0.9,3)/eigenvalues[1:3])
  })
  results[i,] <- colMeans(res)
  # results1[i,] <- colSums(res)
}

round(colMeans(results,na.rm = T),3)
# round(colMeans(results1,na.rm = T),3)
results_5x95xstd <- results
# results_5x9xstd <- results
round(colMeans(results_5x9xstd,na.rm = T),3)
round(colMeans(results_5x95xstd,na.rm = T),3)



results_tmp <- na.omit(results)
colnames(results_tmp) <- c('O','D','PY','PT','L','M','Z')

results_tmp1 <- data.frame(results_tmp[-30,]) %>% pivot_longer(O:Z)

ggplot(results_tmp1)+
  geom_boxplot(aes(y=value, x=name, fill=name),
               position = 'dodge') +
  ylim(c(0,50)) +
  theme_bw() +
  stat_summary(aes(y=value, x=name),
               fun=mean, geom="point", shape=20, size=3, color="red",
               position = position_dodge2(width = 0.75,
                                          preserve = "single")) +
  theme(axis.text = element_text(size=22)) +
  xlab('') + ylab('') +
  guides(fill='none') +
  geom_text(aes(x=name, y=0, label=round(value,3)),
            data=aggregate(.~name, data=results_tmp1, mean),
            col='black', size=3)


# CHeck acf
# Try with oracle detla on D
