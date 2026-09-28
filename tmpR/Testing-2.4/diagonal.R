.fit_diag_model <- function(X,
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
  K_Ce <- which.min(cumsum(eigenvalues)/sum(eigenvalues) < 1-w)

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


    #############
    # METHOD 1
    Chat1 <- Dhat1 <- matrix(0, nrow=dimen, ncol=dimen)
    for(k in 1:(N-1)){
      diag_X2 <- diag_more <- rep(NA, dimen)

      Xk_mx <- X2[,,k] - m_X2
      Z <- ( X2[,,k+1] - m_X2 ) %*% B / r # Tried all transposed of B

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
         'D_hat1'= Dhat1,
         'm_X2'=m_X2
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

  alpha_mat <- matrix(0,n,n)
  alpha_mat[1:nrow(alpha),1:ncol(alpha)] <- alpha

  list('alpha'=diag(alpha),
       # 'delta'=pmax(0,tmp$m_X2 - (diag(diag(alpha),n) %*% tmp$m_X2) / n),
       'delta'=pmax(0,
                    (tmp$m_X2 - (diag(diag(alpha),n) %*% C_eps)) %*% C_eps_dag / n),
       # 'delta1'=pmax(0,
       #              (tmp$m_X2 - (diag(diag(alpha),n) %*% C_eps)) %*% C_eps_dag),
       'delta2'=pmax(0,
                    (tmp$m_X2 - (alpha_mat %*% C_eps)) %*% C_eps_dag / n),
       'K_Ce'=K_Ce)
}

.pred_diag <- function(resolution, alphas, basis, delta, X2_prev, C_eps){
  # Get Alpha
  alpha_applied <- matrix(0, nrow=resolution^2, ncol = resolution^2)
  for(k in 1:length(alphas)){
    hkk <- as.vector(basis[,k] %o% basis[,k])
    alpha_applied <- alpha_applied + alphas[k] * (hkk %o% hkk)
  }

  # Get Sigma
  sigma <- as.vector(delta) + operator(alpha_applied, X2_prev)

  # Get Sqrt(Sigma)
  sigma_svd <- svd(matrix(sigma,
                          ncol=resolution,nrow=resolution))
  sigma_sqrt <- sqrt(resolution) *
    sigma_svd$u %*% diag(sqrt(sigma_svd$d)) %*% t(sigma_svd$v)

  # Predict
  sigma_sqrt %*% C_eps %*% sigma_sqrt
}
