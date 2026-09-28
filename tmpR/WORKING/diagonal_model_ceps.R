### SETUP CODE WITH KNOWN INFO ###
#     -- WORKING CODE --         #
###################
# FUNCTIONS

.getVariables_full <- function(alpha_dimen, resolution, C_eps, scale){

  #   Basis/Eigenvalues based on full C_eps
  eig_Ce <- eigen( C_eps )
  basis <- eig_Ce$vectors * sqrt(resolution)
  eigenvalues <- eig_Ce$values / resolution

  #   Alpha is based on tensor sum of `alpha_dimen` basis functions
  alpha <- matrix(0,nrow=resolution^2,ncol = resolution^2)
  for(k in 1:alpha_dimen){
    hkk <- as.vector(basis[,k] %o% basis[,k])
    alpha <- alpha + scale[k] * (hkk %o% hkk)
  }

  #   Delta is full C_eps
  delta <- C_eps

  #   Return Values
  list('alpha'=alpha,
       'delta'=delta,
       'full'=list('C_eps'=C_eps,'eigenvalues'=eigenvalues,
                   'basis'=basis)
  )
}

generateFARCH <- function(N, resolution, alpha, delta, basis, C_eps,
                          burnin=100, p=1){
  times <- seq(0, 1, length.out=resolution)

  NBurn <- burnin + N
  sim_sigma2_matrix1 <-
    sim_sigma2_matrix2 <-
    matrix(0, ncol = NBurn, nrow = resolution^2)
  arch_data1 <- arch_data2 <-
    matrix(NA, nrow = resolution, ncol = NBurn)

  eps <- t(MASS::mvrnorm(n = NBurn,
                         mu = c(rep(0,resolution)),
                         Sigma = C_eps,
                         empirical = FALSE))

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
    # Sigma2_sqrt1 <- sqrt(resolution) *
    #   Sigma2_svd1$u %*% diag(sqrt(Sigma2_svd1$d)) %*% t(Sigma2_svd1$v)
    Sigma2_sqrt2 <- sqrt(resolution) *
      Sigma2_svd2$u %*% diag(sqrt(Sigma2_svd2$d)) %*% t(Sigma2_svd2$v)
    # arch_data1[,i] <- operator(Sigma2_sqrt1, eps[,i])
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


.estFit <- function(X,
                    C_eps, basis, eigenvalues,
                    w=0.05, p=1, dimen=3){
  #################
  ## Vars
  #################
  N <- ncol(X)
  n <- nrow(X)

  #################
  ## Compute Matrices
  #################
  ## C_eps and related
  K_Ce <- which.min(cumsum(eigenvalues[1:dimen])/sum(eigenvalues[1:dimen]) < 1-w)

  C_eps_dag <- matrix(0, nrow=nrow(C_eps), ncol=ncol(C_eps))
  for(i in 1:K_Ce){
    C_eps_dag <- C_eps_dag +
      ( 1 / eigenvalues[i] ) * (basis[,i] %o% basis[,i] )
  }
  # Note C_eps %*% C_eps_dag for max K_Ce is n^2*Identity

  #################
  ## C/D
  #################
  .CD_hat_new <- function(X, dimen, C_eps_dag, basis){
    N <- ncol(X)
    r <- nrow(X)
    Chat <- Dhat <- matrix(0, nrow=dimen, ncol=dimen)

    #############
    # Common
    B <- C_eps_dag #%*% basis  # Tried all transposed of basis and div by r

    #############
    # METHOD 1
    X2 <- array(apply(X, 2, tcrossprod),
                dim=c(nrow(X),nrow(X),ncol(X)))
    m_X2 <- rowMeans(X2, dims = 2)

    # #############
    # # METHOD 2
    # X_eig <- t( t(X) %*% basis ) / r
    # m_X_eig <- rowMeans(X_eig)
    # X_eig2 <- array(apply(X_eig, 2, tcrossprod),
    #                 dim=c(nrow(X),nrow(X),ncol(X)))
    # m_X_eig2 <- rowMeans(X_eig2, dims = 2)

    #############
    # METHOD 1
    Chat1 <- Dhat1 <- matrix(0, nrow=dimen, ncol=dimen)
    for(k in 1:(N-1)){
      diag_X2 <- diag_more <- rep(NA, dimen)

      Xk_mx <- X2[,,k] - m_X2
      Z <- ( X2[,,k+1] - m_X2 ) %*% B / r

      for(i in 1:dimen){
        diag_X2[i] <- t( Xk_mx %*% basis[,i]) %*% basis[,i] / r^2

        diag_more[i] <- t( Z %*% basis[,i] ) %*% basis[,i] / r^2
      }
      Chat1 <- Chat1 + diag_X2 %*% t(diag_X2)
      Dhat1 <- Dhat1 + diag_X2 %*% t(diag_more)
    }
    # Add N
    diag_X2 <- rep(NA, dimen)
    Xk_mx <- X2[,,N] - m_X2
    for(i in 1:dimen){
      diag_X2[i] <- t( Xk_mx %*% basis[,i]) %*% basis[,i] / r^2
    }
    Chat1 <- Chat1 + diag_X2 %*% t(diag_X2)

    Chat1 <- Chat1 / N
    Dhat1 <- Dhat1 / N

    list('C_hat'= Chat1,
         'D_hat1'= Dhat1
    )
  }

  tmp <- .CD_hat_new(X, dimen, C_eps_dag, basis)

  alpha <- tryCatch({
    tmp$D_hat %*% solve(tmp$C_hat)
  }, error = function(e){
    warning('C_hat not invertible')
    eigs <- eigen(tmp$C_hat)
    dimen_inv <- which.min(cumsum(eigs$values)/sum(eigs$values) < 1-w)
    C_hat_dag <- matrix(0, nrow=nrow(eigs$vectors), ncol = nrow(eigs$vectors))
    for(i in 1:dimen_inv){
      C_hat_dag <- C_hat_dag +
        ( 1 / eigs$values[i] ) * (eigs$vectors[,i] %o% eigs$vectors[,i] )
    }

    tmp$D_hat %*% C_hat_dag
  })

  diag(alpha)
}

run_all_full <- function(N=500, resolution=20,
                         alpha_dimen=2, est_dimen=2,
                    scale=c(0.8,0)){
  ## Key Paramters
  #   Basis
  if(length(scale)!= alpha_dimen){
    if(length(scale)>alpha_dimen) {
      scale <- scale[1:alpha_dimen]
    }else{
      scale <- c(scale, rep(0,alpha_dimen-length(scale)) )
    }
  }
  full_vars <- .getVariables_full(
    alpha_dimen=alpha_dimen,
    resolution = resolution, C_eps = .generate_ceps(resolution),
    scale = scale)


  # Generate data
  data <- generateFARCH(N = N, resolution = resolution,
                        alpha=full_vars$alpha, delta=full_vars$delta,
                        basis=full_vars$full$basis, C_eps=C_eps,
                        burnin=100, p=1)

  # Estimate data
  estimation <- .estFit(X = data$data1,
                        C_eps = full_vars$full$C_eps,
                        basis = full_vars$full$basis,
                        eigenvalues = full_vars$full$eigenvalues,
                        w=0.1, p=1, dimen=est_dimen)

  # Previous project onto the phi
  list('Estimated' = round(estimation, 5),
       'Previous' = round(acf(colMeans(data$data1)^2,plot = FALSE)$acf[1:est_dimen + 1], 5),
       'Data' = data)
}

# Changes
#     Alpha diagional of matrix
#     Eigenvalue fixed

res <- run_all_full(alpha_dimen=2, est_dimen=3,scale = c(0, 0.8))
res$Estimated
res$Previous

## USE OU COV
###################
res <- run_all_full(alpha_dimen=2, est_dimen=3,scale = c(0.8))
res$Estimated
res$Previous

res <- run_all_full(alpha_dimen=2, est_dimen=3,scale = c(0.6, 0.4))
res$Estimated
res$Previous

res <- run_all_full(alpha_dimen=2, est_dimen=20,scale = c(0.6, 0.4))
res$Estimated
res$Previous

res <- run_all_full(alpha_dimen=3, est_dimen=5,scale = c(0, 0, 0.8))
res$Estimated
res$Previous

###################
set.seed(123)
res <- run_all_full(alpha_dimen=2, est_dimen=3,scale = c(0, 0.8))
res_plot <- plot_fd(res$Data$data1)
plotly::save_image(res_plot, 'C:/Users/jerem/Downloads/ceps_d2_e3.png')

nSims <- 1000
sim_res <- sim_res1 <- matrix(NA, nrow=nSims,ncol = 3)

set.seed(123)
for(i in 1:nSims){
  cat(i,', ')
  res <- run_all_full(alpha_dimen=2, est_dimen=3,scale = c(0, 0.8))
  sim_res[i,] <- res$Estimated
  sim_res1[i,] <- res$Previous
}
ggplot() +
  geom_density(aes(x=sim_res[,3],col='New',fill='New'),alpha=0.75 ) +
  geom_density(aes(x=sim_res1[,3],col='Proj',fill='Proj'),alpha=0.75 ) +
  geom_vline(aes(xintercept = 0), linewidth = 1, color='red', linetype = 'dashed') +
  ylab(NULL) +
  xlab(NULL) +
  theme_bw() +
  theme(axis.text = element_text(size=18),
        axis.text.y = element_blank(),
        legend.position = c(0.1,0.8)) +
  scale_color_discrete(name = "") +
  scale_fill_discrete(name = "")
# ceps_d2_e3_term1
