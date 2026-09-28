## FIT MODEL
.operatorARCH <- function(X, w=0.05, thetas=NA, training=0.8, p=1){
  X <- as.matrix(X)
  N <- ncol(X)
  n <- nrow(X)

  # Make tensor square data to use (n x n x N)
  X2_old <- array(apply(X, 2, tcrossprod), dim=c(n,n,N))
  X2 <- matrix(nrow=n^2,ncol=N) # (nn x N)
  for(i in 1:N){
    X2[,i] <- as.vector(X2_old[,,i])
  }

  # C_eps (n x n)
  C_eps <- matrix(NA, n, n)
  times <- 1:n/n
  for (i in 1:n){
    C_eps[i,] <- exp(-abs(times[i]-times)/2)
  }

  #################
  ## USE CV FOR MODEL SETUP
  #     Fit on part, test on pred
  #################
  if(training < 1){
    X_idx_fit <- 1:(training * N)
    # lag back p so can test on the first new data!
    X_idx_pred <- (training * N+1-p):N
  }else{
    X_idx_fit <- 1:N
    X_idx_pred <- 1:N
  }
  N_fit <- length(X_idx_fit)

  # After Eq 23
  #   (n x n)
  C_X_fit <- rowMeans(X2[,X_idx_fit])  # nn

  X2_fit <- X2[,X_idx_fit] # (nn x N_fit)
  X2_pred <- X2[,X_idx_pred] # (nn x N_pred)

  # Define X2_eps (after Eq 17)
  #   (n x n x N_fit), (nn x N_fit)
  X2_eps_fit <- X2_fit
  for(i in 1:N_fit){
    # X2_eps_fit[,i] <- as.vector(adjacent(matrix(X2_fit[,i], ncol=n,nrow=n), C_eps))
    X2_eps_fit[,i] <- as.vector(t(X[, i] %*% C_eps) %*% X[,i])/n^2
  }

  # Define X2p and X2p_eps (Eqs 18)
  #   (pn x n x N_fit), (pnn x N_fit)
  X2p_fit <- X2p_eps_fit <- matrix(nrow=p*n*n, ncol=N_fit)
  for(k in 1:N_fit){
    tmp <- tmp_eps <- list()
    for(j in 1:p){
      if(k-j+1 >0){ # We can look back at the lag
        tmp[[j]] <- X2_fit[,k-j+1,drop=FALSE]
        tmp_eps[[j]] <- X2_eps_fit[,k-j+1,drop=FALSE]
      } else{ # Too early in the sample to look back to data
        tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n^2, ncol=1)
      }
    }
    X2p_fit[,k] <- do.call("rbind", tmp)
    X2p_eps_fit[,k] <- do.call("rbind", tmp_eps)
  }

  # After Eq 23
  #   (pn x n)
  mxp <- matrix(0, ncol=1, nrow=nrow(X2p_fit))
  for(l in p:N_fit){
    mxp <- mxp + X2p_fit[,l]
  }
  mxp <- mxp / (N_fit-p+1)

  # Eqs 20 and 21
  #   C: (pnn x pnn),
  #   D: (nn x pnn)
  tmp <- .CD_hat(X = X2_fit, Xp=X2p_fit, Xpe = X2p_eps_fit, p = p)
  C_hat <- tmp$C_hat
  D_hat <- tmp$D_hat

  # Around Eq 22-23
  #   projCC: (pnn x pnn)
  #   projCe: (n x n)
  tmp <- .computeProjectorsAndInitialize(C_hat, C_eps, w, n,p)
  projector_CC <- tmp$projector_CC # rows of C_hat
  projector_Ce <- tmp$projector_Ce # res of X2
  guess_theta <- tmp$theta
  guess_vu <- tmp$vu


  # Find the optimal Tikhonov parameters
  if(length(thetas)!=2 ||
     sum(is.na(thetas)) > 0 ){

    tryCatch({
      optim_results <-
        optim(c(guess_theta,guess_vu), fn=function(theta0, C_hat, D_hat,
                                                   projector_CC, projector_Ce,
                                                   n, p1, X21, C_X1, mxp){
          theta <- theta0[1]
          vu <- theta0[2]

          ests <- .estim_params_ARCH_op(theta = theta, vu=vu, C_hat = C_hat,
                                        D_hat = D_hat, C_eps=C_eps, C_X=C_X1,
                                        projector_CC = projector_CC,
                                        projector_Ce = projector_Ce,
                                        n = n, p = p1, mxp = mxp)

          .compute_err_ARCH_op(
            X2 = X21, alpha = ests$alpha_hat, delta = ests$delta_hat,
            C_eps=C_eps, p=p1, C_X=C_X1)
        },
        method="L-BFGS-B",
        lower=c(10^-8,0), upper=c(Inf,Inf),
        C_hat=C_hat, D_hat=D_hat, projector_CC=projector_CC,
        projector_Ce=projector_Ce, n=n, p1=p, X21=X2_pred,
        C_X1=C_X_fit, mxp = mxp)
    }, error =function(e){
      optim_results <-
        optim(c(guess_theta*2,guess_vu*2), fn=function(theta0, C_hat, D_hat,
                                                       projector_CC, projector_Ce,
                                                       n, p1, X21, C_X1, mxp){
          theta <- theta0[1]
          vu <- theta0[2]

          ests <- .estim_params_ARCH_op(theta = theta, vu=vu, C_hat = C_hat,
                                        D_hat = D_hat, C_eps=C_eps, C_X=C_X1,
                                        projector_CC = projector_CC,
                                        projector_Ce = projector_Ce,
                                        n = n, p = p1, mxp=mxp)

          .compute_err_ARCH_op(
            X2 = X21, alpha = ests$alpha_hat, delta = ests$delta_hat,
            C_eps=C_eps, p=p1, C_X=C_X1)
        },
        method="L-BFGS-B",
        lower=c(10^-8,0),upper=c(Inf,Inf),
        C_hat=C_hat, D_hat=D_hat, projector_CC=projector_CC,
        projector_Ce=projector_Ce, n=n, p1=p, X21=X2_pred,
        C_X1=C_X_fit, mxp = mxp)


    })

    thetas <- optim_results$par
  }

  #################
  ## USE ALL DATA FOR FINAL SETUP
  #################

  C_X <- rowMeans(X2)#,dims=2)
  # Define X2_eps
  X2_eps <- X2
  for(i in 1:N){
    X2_eps[,i] <- as.vector(adjacent(matrix(X2[,i], ncol=n,nrow=n), C_eps))
    # X2_eps[,i] <- as.vector(matrix(X2[,i], ncol=n,nrow=n) %*% C_eps) / n
  }

  # Define X2p
  X2p <- X2p_eps <- matrix(nrow=p*n*n, ncol=N)
  for(k in 1:N){
    tmp <- tmp_eps <- list()
    for(j in 1:p){
      if(k-j+1 >0){ # Lag back if enough data
        tmp[[j]] <- X2[,k-j+1,drop=FALSE]
        tmp_eps[[j]] <- X2_eps[,k-j+1,drop=FALSE]
      } else{ # Can't lag back (too early)
        tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n^2, ncol=1)
      }
    }
    X2p[,k] <- do.call("rbind", tmp)
    X2p_eps[,k] <- do.call("rbind", tmp_eps)
  }

  mxp <- matrix(0, ncol=1, nrow=nrow(X2p))
  for(l in p:N){
    mxp <- mxp + X2p[,l]
  }
  mxp <- mxp / (N-p+1)

  # Eqs 35 and 36
  tmp <- .CD_hat(X = X2, Xp=X2p, Xpe = X2p_eps, p = p)
  C_hat <- tmp$C_hat
  D_hat <- tmp$D_hat

  tmp <- .computeProjectorsAndInitialize(C_hat, C_eps,w,n,p)
  projector_CC <- tmp$projector_CC
  projector_Ce <- tmp$projector_Ce

  ests <- .estim_params_ARCH_op(theta = thetas[1], vu=thetas[2],
                                C_hat = C_hat,
                                D_hat = D_hat, C_eps=C_eps, C_X=C_X,
                                projector_CC = projector_CC,
                                projector_Ce = projector_Ce,
                                n = n, p = p, mxp = mxp)

  # .compute_err_ARCH_op(
  #   X2 = X2, alpha = alpha_true_tmp, delta = (sapply(times,delta_par,t=times)),
  #   C_eps=C_eps, p=p, C_X=C_X_fit)
  # .compute_err_ARCH_op(
  #   X2 = X2, alpha = ests$alpha_hat, delta = ests$delta_hat,
  #   C_eps=C_eps, p=p, C_X=C_X)


  list('hat_alpha'=ests$alpha_hat,
       'hat_delta'=ests$delta_hat,
       'fitting_params'=list('thetas'=thetas))
}

################


.computeProjectorsAndInitialize <- function(C_hat, C_eps,w,n,p){
  # Between Eqs 22-23
  eig_CC <- eigen( adjacent(C_hat,t(Conj(C_hat)) )) # C for p=1
  eig_Ce <- eigen( C_eps )
  CC_n <- nrow(eig_CC$vectors)
  Ce_n <- nrow(eig_Ce$vectors)
  eig_CC$vectors <- eig_CC$vectors * sqrt(CC_n)
  eig_Ce$vectors <- eig_Ce$vectors * sqrt(Ce_n)

  K_CC <- which.min(cumsum(eig_CC$values)/sum(eig_CC$values) < 1-w)
  K_Ce <- which.min(cumsum(eig_Ce$values)/sum(eig_Ce$values) < 1-w)

  # Compute Projection (TODO:: Verify)
  #   CC: (pnn x pnn) = (pnn x pnn)x(pnn x pnn)
  #   Ce: (n x n)
  projector_CC <- matrix(0, nrow=CC_n, ncol=CC_n)
  projector_Ce <- matrix(0, nrow=Ce_n, ncol=Ce_n)

  for(i in 1:K_CC){
    # (pnn x pnn)
    projector_CC <- projector_CC +
      tensor_vec( eig_CC$vectors[,i], eig_CC$vectors[,i] )
  }
  for(i in 1:K_Ce){
    # (n x n)
    projector_Ce <- projector_Ce +
      tensor_vec( eig_Ce$vectors[,i], eig_Ce$vectors[,i] )
  }

  ## Initialize
  theta <- eig_CC$value[1] / 2
  vu <- eig_Ce$value[1] / 2

  list(projector_CC=projector_CC,
       projector_Ce=projector_Ce,
       theta=theta,
       vu=vu)
}

.estim_params_ARCH_op <- function(theta, vu, C_hat,D_hat, C_eps, C_X,
                                  projector_CC,projector_Ce, n, p, mxp){
  # Compute alpha, Eq 22
  C_hat_st <- t(Conj(C_hat)) # (pnn x pnn)
  C_dag <-
    adjacent(C_hat_st,
             solve( adjacent(C_hat, C_hat_st) +
                      diag(theta,nrow(C_hat)) ) # (pnn x pnn)
    )
  # (nn x pnn) x (pnn x pnn) x (pnn x pnn)
  alpha_hat <- adjacent(adjacent(D_hat, C_dag), projector_CC)
  #alpha_hat  <- apply(alpha_hat, 2, pmax, 0) # nn x nnp

  ## Delta, Eq 39
  C_eps_diag <- solve( C_eps + diag(vu, n) ) # n x n
  alpha_sum <- matrix(operator(alpha_hat, mxp),
                      nrow = sqrt(nrow(alpha_hat)),
                      ncol = sqrt(nrow(alpha_hat)) )

  delta_hat <-
    adjacent(
      adjacent(
        C_X - adjacent(alpha_sum, C_eps),
        C_eps_diag ),
      projector_Ce)

  delta_hat <- apply(delta_hat, 2, pmax, 0)

  list('alpha_hat'=alpha_hat,
       'delta_hat'=delta_hat)
}

##################


.compute_err_ARCH_op <- function(X2, alpha, delta, C_eps, C_X, p){
  error <- rep(0, ncol(X2))
  estimate <- .pred_ARCH_op(X2, alpha, delta, C_eps, p)

  for(i in (1+p):ncol(X2) ){
    error[i] <- dot_integrate((X2[,i] - estimate[,i])^2)
  }

  sum(error)
}

.pred_ARCH_op <- function(X2, alpha, delta, C_eps, p){
  N <- ncol(X2)
  n <- nrow(delta)
  Sigma2 <- array(NA,dim=c(n,n,N))
  Estimate <- matrix(NA, nrow = n^2,ncol = N)

  # Eq 2
  for(i in (1+p):N){
    tmp_op <- matrix(0, nrow = n^2,ncol = 1)
    for(j in 1:p){
      ## TODO: Question
      tmp_op <- tmp_op + operator(alpha,X2[,i-j])
      # matrix(alpha[,1:nrow(alpha)+(j-1)*nrow(alpha)] %*% X2[,i-j],
      #        ncol=n, nrow=n)
    }
    Sigma2[,,i] <- delta + matrix(tmp_op,nrow = n,ncol = n)

    Sigma_svd <- tryCatch({ svd(matrix(Sigma2[,,i],n,n)) },
                          error=function(e){ svd::propack.svd(matrix(Sigma2[,,i],n,n)) })
    Sigma_sqrt <- sqrt(n) * Sigma_svd$u %*%
      diag(sqrt(Sigma_svd$d)) %*%
      t(Sigma_svd$v)

    ## Eq 3
    Estimate[,i] <- adjacent(adjacent(Sigma_sqrt, C_eps), Sigma_sqrt) #/ n^2
  }

  Estimate
}
