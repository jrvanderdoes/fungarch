generateFARCH_use <- function(N, resolution, alpha, delta, basis,
                              burnin=100, p=1, eigs=NULL,
                              C_eps=NULL){
  #### Functions
  .get_kl_errors <- function(N, resolution, basis, eigs){

    .get_errors <- function(basis,eigs){
      ## NEXT
      xi <- sapply(eigs,
                   function(e) {
                     stats::rnorm(1, mean = 0, sd = 1)*sqrt(e)
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

  if(is.null(eigs)){
    eps <- t(MASS::mvrnorm(n = NBurn,
                           mu = c(rep(0,resolution)),
                           Sigma = C_eps,
                           empirical = FALSE))
  }else{
    eps <- .get_kl_errors(N = NBurn, resolution = resolution,
                          basis=basis, eigs=eigs)
  }

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
