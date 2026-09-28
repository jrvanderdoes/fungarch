# Moore Penrose
.est_opARCH_mp <- function(X, w=0.05, p=1){
  N <- ncol(X)
  n <- nrow(X)

  # Make tensor square data to use (n x n x N)
  X2_old <- array(apply(X, 2, tcrossprod), dim=c(nrow(X),nrow(X),ncol(X)))
  X2 <- matrix(nrow=n^2,ncol=N) # (nn x N)
  for(i in 1:N){
    X2[,i] <- as.vector(X2_old[,,i])
  }

  # C_eps (n x n)
  C_eps <- .generate_ceps(n)

  eig_Ce <- eigen( C_eps )
  Ce_n <- nrow(eig_Ce$vectors)
  eig_Ce$vectors <- eig_Ce$vectors * sqrt(Ce_n)
  eig_Ce$values <- eig_Ce$values / Ce_n # ADDED

  K_Ce <- which.min(cumsum(eig_Ce$values)/sum(eig_Ce$values) < 1-w)

  C_eps_dag <- matrix(0, nrow=nrow(C_eps), ncol=ncol(C_eps))
  for(i in 1:K_Ce){
    C_eps_dag <- C_eps_dag +
      ( 1 / eig_Ce$values[i] ) * (eig_Ce$vectors[,i] %o% eig_Ce$vectors[,i] )
  }

  #################
  ## Fit Model
  #################
  C_X <- as.vector(adjacent(matrix(rowMeans(X2),n,n), C_eps))

  # Define X2p
  X2p <- matrix(nrow=p*n*n, ncol=N)
  for(k in 1:N){
    tmp <- list()
    for(j in 1:p){
      if(k-j+1 >0){ # Lag back if enough data
        tmp[[j]] <- X2[,k-j+1,drop=FALSE]
        # tmp_eps[[j]] <- X2_eps[,k-j+1,drop=FALSE]
      } else{ # Can't lag back (too early)
        tmp[[j]] <- matrix(0, nrow=n^2, ncol=1)
      }
    }
    X2p[,k] <- do.call("rbind", tmp)
    # X2p_eps[,k] <- do.call("rbind", tmp_eps)
  }

  # Eqs 35 and 36
  tmp <- .CD_hat(X2 = X2, X2p=X2p, p = p, C_eps_dag=C_eps_dag)
  C_hat <- tmp$C_hat
  D_hat <- tmp$D_hat

  # # After Eq 36 and Eq 39
  eig_C <- eigen( adjacent(C_hat,t(Conj(C_hat)) ) )
  CC_n <- nrow(eig_C$vectors)
  # Make Correct
  eig_C$vectors <- eig_C$vectors * sqrt(CC_n)
  eig_C$values <- eig_C$values / CC_n # ADDED

  # tmp=fda::pca.fd(fda::Data2fd( adjacent(C_hat,t(Conj(C_hat)) ) ),
  #                 nharm = ncol(C_hat)-1)
  # eig_C$vectors <- tmp$scores
  # eig_C$values <- tmp$values
  # tmp=fda::pca.fd(fda::Data2fd( C_eps ),
  #                 nharm = n)
  # eig_Ce$vectors <- tmp$scores
  # eig_Ce$values <- tmp$values

  # pca_C <- pca_func( adjacent( C_hat,t(Conj(C_hat)) ), TVE=1-w )
  # eig_C <- list('values'=pca_C$sdev^2,
  #               'vectors'=pca_C$rotation )
  # pca_Ce <- pca_func( C_eps, TVE=1-w )
  # eig_Ce <- list('values'=pca_Ce$sdev^2,
  #                'vectors'=pca_Ce$rotation )
  # K_C <- length(eig_C$values)

  K_C <- which.min(cumsum(eig_C$values)/sum(eig_C$values) < 1-w)

  # Estimate alpha
  C_dag <- matrix(0, nrow=nrow(C_hat), ncol=ncol(C_hat))
  for(i in 1:K_C){
    C_dag <- C_dag +
      ( 1 / eig_C$values[i] ) * (eig_C$vectors[,i] %o% eig_C$vectors[,i] )
  }
  # n (np) x np (np)
  alpha_hat <- adjacent(D_hat, C_dag)
  # adjacent does not matter in scaling
  #### max(dat$params$alpha) / max(alpha_hat)

  # alpha_hat <- dat$params$alpha

  ## Delta, Eq 39
  mxp <- .mean(X = X2p,p = p,lag = 0,lag2 = 0)


  alpha_applied <- operator(alpha_hat, mxp)


  delta_hat <-
    apply(adjacent(matrix(C_X,n,n) - adjacent(matrix(alpha_applied,n,n), C_eps), C_eps_dag), 2, pmax, 0)
  # plot_fd(delta_hat)
  #### max(dat$params$delta)/max(delta_hat)


  list('hat_alpha'=alpha_hat,
       'hat_delta'=delta_hat,
       'C_hat'=C_hat,
       'D_hat'=D_hat)
}

.CD_hat <- function(X2, X2p, p, C_eps_dag){

  N <- ncol(X2p)

  # Eqs 20 and 21
  #   C: (pnn x pnn)
  #   D: (nn x pnn)
  Chat <-  matrix(0, nrow=nrow(X2p), ncol=nrow(X2p))
  Dhat <-  matrix(0, nrow=nrow(X2), ncol=nrow(X2p))

  # After Eq 35
  mx2p <- .mean(X2p,p=p,lag=0,lag2=0)
  mx2p_pr <- .mean(X2p,p=p,lag=1,lag2=0)
  mx2_pr <- .mean(X2,p=p,lag=1,lag2=1)

  Chat <- tensor_vec(X2[,N]-mx2p, X2p[,N]-mx2p)
  for(k in p:(N-1)){
    ## C_hat, Eq 20
    #   N case computed before loop
    Chat <- Chat + tensor_vec(X2[,k]-mx2p, X2p[,k]-mx2p)

    ## D_hat, Eq 21
    Dhat <- Dhat + tensor_vec(X2p[,k]-mx2p_pr,
                              as.vector(
                                adjacent(matrix(X2[,k+1]-mx2_pr, sqrt(nrow(X2)), sqrt(nrow(X2))),
                                              C_eps_dag))) # mattric unstacked
  }

  list('C_hat'= Chat / (N-p),
       'D_hat'= Dhat / (N-p-1)
  )
}

.mean <- function(X, p, lag,lag2=0){
  mat <- matrix(0, nrow=nrow(X), ncol=1)
  N <- ncol(X)

  for(k in p:(N-lag)){
    mat <- mat + X[,k+lag2]
  }

  mat / (N-p+1-lag)
}
