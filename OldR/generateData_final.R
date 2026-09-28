
generateFunctionalARCH_final <- function(resolution, N,
                                        delta=.generate_ceps(resolution),
                                        burnin=100, p=1,
                                        K=2, scale=0.4){
  times <- seq(0, 1, length.out=resolution)

  NBurn <- burnin + N
  sim_sigma2_matrix1 <-
  sim_sigma2_matrix2 <-
    matrix(0, ncol = NBurn, nrow = resolution^2)
  arch_data <- arch_data1 <-
  arch_data2 <-
    matrix(NA, nrow = resolution, ncol = NBurn)

  eps <- .generateOU(resolution = resolution, N = NBurn)

  ## Define variance and GARCH matrices
  # Initialize
  for(i in 1:p){
    arch_data1[,i] <- operator(delta, eps[,i])
    arch_data2[,i] <- operator(delta, eps[,i])
  }

  # alpha <- readRDS(alpha_par)

  # Set matrices
  Sigma2_list1 <- Sigma2_sqrt_list1 <-
    Sigma2_list2 <- Sigma2_sqrt_list2 <- list()

  # SETUP BASIS
  val_seq <- seq(0, 1, length.out=resolution)
  basis_vals <- .generate_basis(resolution)

  # ALPHA
  alpha <- matrix(0,nrow=resolution^2,ncol = resolution^2)
  for(k in 1:K){
    hkk <- tensor_vec(basis_vals[,k], basis_vals[,k])
    alpha <- alpha + tensor_vec(hkk, hkk)
  }
  alpha <- scale * alpha


  # Compute values
  for(i in (p+1):NBurn){
    tmp_alpha1 <- tmp_alpha2 <- rep(0, resolution^2)
    for(pp in 1:p){
      # Basis
      prev_tensor <- tensor_vec(arch_data1[,i-pp], arch_data1[,i-pp])
      for(b in 2:K){
        mat_val <- basis_vals[,b] %o% basis_vals[,b]
        tmp <- mat_val * scale *
          as.numeric((basis_vals[,b] %*% arch_data1[,i-pp] / nrow(basis_vals))^2)
          # as.numeric((basis_vals[,b] %*% arch_data1[,i-pp])^2 / nrow(basis_vals)^2)
        tmp_alpha1 <- tmp_alpha1 + tmp
      }

      prev_tensor <- tensor_vec(arch_data2[,i-pp], arch_data2[,i-pp])
      # Alpha
      tmp_alpha2 <- tmp_alpha2 + operator(alpha, prev_tensor)
    }

    #  Sigma^2
    sim_sigma2_matrix1[,i] <-
      as.vector(delta) +
      tmp_alpha1
    sim_sigma2_matrix2[,i] <-
      as.vector(delta) +
      tmp_alpha2

    # Attempt 2: Compute the sqrt matrix
    tryCatch({
      Sigma2_svd1 <- svd(matrix(sim_sigma2_matrix1[,i],
                                ncol=resolution,nrow=resolution))
      Sigma2_svd2 <- svd(matrix(sim_sigma2_matrix2[,i],
                                ncol=resolution,nrow=resolution))
    }, error=function(e){
      cat(i)
      # View(sim_sigma2_matrix1)
      # View(arch_data1)
      stop('Error')
    })
    Sigma2_sqrt1 <- sqrt(resolution) *
      Sigma2_svd1$u %*% diag(sqrt(Sigma2_svd1$d)) %*% t(Sigma2_svd1$v)
    Sigma2_sqrt2 <- sqrt(resolution) *
      Sigma2_svd2$u %*% diag(sqrt(Sigma2_svd2$d)) %*% t(Sigma2_svd2$v)
    arch_data1[,i] <- operator(Sigma2_sqrt1, eps[,i])
    arch_data2[,i] <- operator(Sigma2_sqrt2, eps[,i])

    Sigma2_sqrt_list1[[i]] <- Sigma2_sqrt1
    Sigma2_sqrt_list2[[i]] <- Sigma2_sqrt2
  }

  # Return
  new_res1 <- list(data = arch_data1[, burnin + 1:N],
    data1 = arch_data2[, burnin + 1:N],
    sigma_mat = sim_sigma2_matrix1[, burnin + 1:N],
    sigma_mat1 = sim_sigma2_matrix2[, burnin + 1:N],
    error_matrix = eps[, burnin + 1:N],
    params=list(p=p,q=q,
                # alpha_par=alpha_par,
                alpha = alpha,
                delta = delta)
  )
}


.generateOU <- function(resolution, N){

  times <- 1:resolution/resolution

  # Covariance structure (OU process Cov)
  comat <- matrix(NA,resolution,resolution)
  for (i in 1:resolution){
    comat[i,] <- exp(-abs(times[i]-times)/2)
  }

  fiid <- MASS::mvrnorm(n = N,
                        mu = c(rep(0,resolution)),
                        Sigma = comat,
                        empirical = FALSE)

  t(fiid)
}
