## GENERATE DATA
generateFunctionalARCH <- function(resolution, N,
                                   alpha_par=NULL, delta_par=NULL,
                                   burnin=100, p=1){
  times <- seq(0, 1, length.out=resolution)

  NBurn <- burnin + N
  sim_sigma2_matrix1 <- matrix(0, ncol = NBurn, nrow = resolution^2)
  arch_data <- arch_data1 <-  matrix(NA, nrow = resolution, ncol = NBurn)

  eps <- .generateOU(resolution = resolution, N = NBurn)

  ## Define variance and GARCH matrices
  # Initialize
  for(i in 1:p){
    arch_data1[,i] <- operator(sapply(times,delta_par, t=times), eps[,i])
  }

  alpha <- readRDS(alpha_par)

  # Set matrices
  Sigma2_list1 <- Sigma2_sqrt_list1 <- list()

  # Compute values
  for(i in (p+1):NBurn){
    tmp_alpha1 <- rep(0, resolution^2)
    for(pp in 1:p){
      prev_tensor <- tensor_vec(arch_data1[,i-pp],arch_data1[,i-pp])

      tmp_alpha1 <- tmp_alpha1 + operator(alpha, prev_tensor)
    }

    #  Sigma^2
    sim_sigma2_matrix1[,i] <-
      as.vector(sapply(times,delta_par,t=times)) +
      tmp_alpha1

    # Attempt 2: Compute the sqrt matrix
    tryCatch({
      Sigma2_svd1 <- svd(matrix(sim_sigma2_matrix1[,i],
                                ncol=resolution,nrow=resolution))
    }, error=function(e){
      cat(i)
      View(sim_sigma2_matrix1)
      View(arch_data1)
      stop('Error')
    })
    Sigma2_sqrt1 <- sqrt(resolution) *
      Sigma2_svd1$u %*% diag(sqrt(Sigma2_svd1$d)) %*% t(Sigma2_svd1$v)
    arch_data1[,i] <- operator(Sigma2_sqrt1, eps[,i])

    Sigma2_sqrt_list1[[i]] <- Sigma2_sqrt1
  }

  # Return
  new_res1 <- list(data = arch_data1[, burnin + 1:N],
                   sigma_mat = sim_sigma2_matrix1[, burnin + 1:N],
                   error_matrix = eps[, burnin + 1:N],
                   params=list(p=p,q=q,
                               alpha_par=alpha_par,
                               alpha = alpha,
                               delta_par=delta_par,
                               delta = as.vector(sapply(times,delta_par,t=times)),
                               Sigma2_list= Sigma2_list1[burnin + 1:N],
                               Sigma2_sqrt_list= Sigma2_sqrt_list1[burnin + 1:N])
  )
}


generateFunctionalARCH_basis <- function(resolution, N,
                                         delta_par=NULL,
                                         burnin=100, p=1,
                                         K=4){
  times <- seq(0, 1, length.out=resolution)

  NBurn <- burnin + N
  sim_sigma2_matrix1 <- matrix(0, ncol = NBurn, nrow = resolution^2)
  arch_data <- arch_data1 <-  matrix(NA, nrow = resolution, ncol = NBurn)

  eps <- .generateOU(resolution = resolution, N = NBurn)

  ## Define variance and GARCH matrices
  # Initialize
  for(i in 1:p){
    arch_data1[,i] <- operator(sapply(times,delta_par, t=times), eps[,i])
  }

  # alpha <- readRDS(alpha_par)

  # Set matrices
  Sigma2_list1 <- Sigma2_sqrt_list1 <- list()

  # SETUP BASIS
  val_seq <- seq(0, 1, length.out=resolution)
  tmp <- fda::create.fourier.basis(nbasis = K, rangeval = c(0,1))
  basis_vals <- fda::eval.basis(seq(0,1,length.out=length(val_seq)), tmp)


  # Compute values
  for(i in (p+1):NBurn){
    tmp_alpha1 <- rep(0, resolution^2)
    for(pp in 1:p){

      prev_tensor <- arch_data1[,i-pp] %*% t(arch_data1[,i-pp])
      for(b in 2:K){
        mat_val <- basis_vals[,b] %o% basis_vals[,b]
        tmp <- 20*mat_val *
          as.numeric((basis_vals[,b] %*% arch_data1[,i-pp])^2 / nrow(basis_vals)^2)
        tmp_alpha1 <- tmp_alpha1 + tmp
      }
    }

    #  Sigma^2
    sim_sigma2_matrix1[,i] <-
      as.vector(sapply(times,delta_par,t=times)) +
      tmp_alpha1

    # Attempt 2: Compute the sqrt matrix
    tryCatch({
      Sigma2_svd1 <- svd(matrix(sim_sigma2_matrix1[,i],
                                ncol=resolution,nrow=resolution))
    }, error=function(e){
      cat(i)
      View(sim_sigma2_matrix1)
      View(arch_data1)
      stop('Error')
    })
    Sigma2_sqrt1 <- sqrt(resolution) *
      Sigma2_svd1$u %*% diag(sqrt(Sigma2_svd1$d)) %*% t(Sigma2_svd1$v)
    arch_data1[,i] <- operator(Sigma2_sqrt1, eps[,i])

    Sigma2_sqrt_list1[[i]] <- Sigma2_sqrt1
  }

  # Return
  new_res1 <- list(data = arch_data1[, burnin + 1:N],
                   sigma_mat = sim_sigma2_matrix1[, burnin + 1:N],
                   error_matrix = eps[, burnin + 1:N],
                   params=list(p=p,q=q,
                               # alpha_par=alpha_par,
                               # alpha = alpha,
                               delta_par=delta_par,
                               delta = as.vector(sapply(times,delta_par,t=times)),
                               Sigma2_list= Sigma2_list1[burnin + 1:N],
                               Sigma2_sqrt_list= Sigma2_sqrt_list1[burnin + 1:N])
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
