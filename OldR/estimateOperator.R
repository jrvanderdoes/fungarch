#' Operator ARCH Estimation
#'
#' Estimates the parameters for ARCH(1) model using operator Yule Walker
#'  estimation.
#'
#' @inheritParams pointwiseARCH
#'
#' @return
#' @export
#'
#' @examples
operatorARCH <- function(X, w=0.05, thetas=NA, training=0.8, p=1){
  N <- ncol(X)
  n <- nrow(X)

  # Make tensor square data to use (n x n x N)
  X2_old <- array(apply(X, 2, tcrossprod),
              dim=c(nrow(X),nrow(X),ncol(X)))
  X2 <- matrix(nrow=n^2,ncol=N) # (nn x N)
  for(i in 1:N){
    X2[,i] <- as.vector(X2_old[,,i])
  }

  # C_eps (n x n)
  C_eps <- matrix(NA, n, n)
  times <- 1:n/n
  for (i in 1:n){
    #C_eps[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
    C_eps[i,] <- exp(-abs(times[i]-times)/2)
  }
  # C_eps <- as.vector(C_eps) #(nn)

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
  #C_X_fit <- rowMeans(X2[,,X_idx_fit],dims=2)
  C_X_fit <- rowMeans(X2[,X_idx_fit])  # nn

  # X2_fit <- X2[,,X_idx_fit] # (n x n x N_fit)
  # X2_pred <- X2[,,X_idx_pred] # (n x n x N_pred)
  X2_fit <- X2[,X_idx_fit] # (nn x N_fit)
  X2_pred <- X2[,X_idx_pred] # (nn x N_pred)

  # Define X2_eps (after Eq 17)
  #   (n x n x N_fit), (nn x N_fit)
  X2_eps_fit <- X2_fit
  for(i in 1:N_fit){
    # X2_eps_fit[,,i] <- (X2_fit[,,i] %*% C_eps) / n
    X2_eps_fit[,i] <- as.vector(matrix(X2_fit[,i], ncol=n,nrow=n) %*% C_eps) / n
  }

  # Define X2p and X2p_eps (Eqs 18)
  #   (pn x n x N_fit), (pnn x N_fit)
  # X2p_fit <- X2p_eps_fit <- array(dim=c(n*p,n,N))
  X2p_fit <- X2p_eps_fit <- matrix(nrow=p*n*n, ncol=N_fit)
  for(k in 1:N_fit){
    tmp <- tmp_eps <- list()
    for(j in 1:p){
      if(k-j+1 >0){ # We can look back at the lag
        # tmp[[j]] <- X2_fit[,,k-j+1]
        # tmp_eps[[j]] <- X2_eps_fit[,,k-j+1]
        tmp[[j]] <- X2_fit[,k-j+1,drop=FALSE]
        tmp_eps[[j]] <- X2_eps_fit[,k-j+1,drop=FALSE]
      } else{ # Too early in the sample to look back to data
        # tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n, ncol=n)
        tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n^2, ncol=1)
      }
    }
    # X2p_fit[,,k] <- t(do.call("cbind", tmp))
    # X2p_eps_fit[,,k] <- t(do.call("cbind", tmp_eps))
    X2p_fit[,k] <- do.call("rbind", tmp)
    X2p_eps_fit[,k] <- do.call("rbind", tmp_eps)
  }

  # After Eq 23
  #   (pn x n)
  # mxp <- as.matrix(array(0,dim(X2p_fit[,,1])))
  mxp <- matrix(0, ncol=1, nrow=nrow(X2p_fit))
  for(l in p:N_fit){
    # mxp <- mxp + X2p_fit[,,l]
    mxp <- mxp + X2p_fit[,l]
  }
  mxp <- mxp / (N_fit-p+1)

  # Eqs 20 and 21
  #   C: (pnn x pnn),
  #   D: (nn x pnn)
  tmp <- .estim_CD_hat(X = X2_fit, Xp=X2p_fit, Xpe = X2p_eps_fit, p = p)
  C_hat <- tmp$C_hat
  D_hat <- tmp$D_hat

  # Around Eq 22-23
  #   projCC: (pnn x pnn)
  #   projCe: (n x n)
  tmp <- .computeProjectorsAndInitialization(C_hat, C_eps, w, n,p)
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

          ests <- .estimate_params_ARCH_op(theta = theta, vu=vu, C_hat = C_hat,
                                           D_hat = D_hat, C_eps=C_eps, C_X=C_X1,
                                           projector_CC = projector_CC,
                                           projector_Ce = projector_Ce,
                                           n = n, p = p1, mxp = mxp)

          .compute_ARCH_error_op(
            X2 = X21, alpha = ests$alpha_hat, delta = ests$delta_hat,
            C_eps=C_eps, p=p1, C_X=C_X1)
        },
        method="L-BFGS-B",#"Brent",
        lower=c(10^-8,0),upper=c(Inf,Inf),
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

          ests <- .estimate_params_ARCH_op(theta = theta, vu=vu, C_hat = C_hat,
                                           D_hat = D_hat, C_eps=C_eps, C_X=C_X1,
                                           projector_CC = projector_CC,
                                           projector_Ce = projector_Ce,
                                           n = n, p = p1, mxp=mxp)

          .compute_ARCH_error_op(
            X2 = X21, alpha = ests$alpha_hat, delta = ests$delta_hat,
            C_eps=C_eps, p=p1, C_X=C_X1)
        },
        method="L-BFGS-B",#"Brent",
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

  # N <- ncol(X) # Number of observations
  C_X <- rowMeans(X2)#,dims=2)
  # Define X2_eps
  X2_eps <- X2
  for(i in 1:N){
    # X2_eps[,,i] <- X2[,,i] %*% C_eps / nrow(X2[,,i])
    X2_eps[,i] <- as.vector(matrix(X2[,i], ncol=n,nrow=n) %*% C_eps) / n
  }

  # Define X2p
  #X2p <- X2p_eps <- array(dim=c(n*p,n,N))
  X2p <- X2p_eps <- matrix(nrow=p*n*n, ncol=N)
  for(k in 1:N){
    tmp <- tmp_eps <- list()
    for(j in 1:p){
      if(k-j+1 >0){ # Lag back if enough data
        # tmp[[j]] <- X2[,,k-j+1]
        # tmp_eps[[j]] <- X2_eps[,,k-j+1]
        tmp[[j]] <- X2[,k-j+1,drop=FALSE]
        tmp_eps[[j]] <- X2_eps[,k-j+1,drop=FALSE]
      } else{ # Can't lag back (too early)
        # tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n, ncol=n)
        tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n^2, ncol=1)
      }
    }
    # X2p[,,k] <- t(do.call("cbind", tmp))
    # X2p_eps[,,k] <- t(do.call("cbind", tmp_eps))
    X2p[,k] <- do.call("rbind", tmp)
    X2p_eps[,k] <- do.call("rbind", tmp_eps)
  }

  # mxp <- matrix(array(0,dim(X2p[,,1])))
  mxp <- matrix(0, ncol=1, nrow=nrow(X2p))
  for(l in p:N){
    # mxp <- mxp + X2p[,,l]
    mxp <- mxp + X2p[,l]
  }
  mxp <- mxp / (N-p+1)

  # Eqs 35 and 36
  tmp <- .estim_CD_hat(X = X2, Xp=X2p, Xpe = X2p_eps, p = p)
  C_hat <- tmp$C_hat
  D_hat <- tmp$D_hat

  tmp <- .computeProjectorsAndInitialization(C_hat, C_eps,w,n,p)
  projector_CC <- tmp$projector_CC
  projector_Ce <- tmp$projector_Ce

  ests <- .estimate_params_ARCH_op(theta = thetas[1], vu=thetas[2],
                                   C_hat = C_hat,
                                   D_hat = D_hat, C_eps=C_eps, C_X=C_X,
                                   projector_CC = projector_CC,
                                   projector_Ce = projector_Ce,
                                   n = n, p = p, mxp = mxp)

  list('hat_alpha'=ests$alpha_hat,
       'hat_delta'=ests$delta_hat,
       'error'=.compute_ARCH_error_op(
         X2 = X2, alpha = ests$alpha_hat, delta = ests$delta_hat,
         C_eps=C_eps, p=p, C_X=C_X),
       'fitting_params'=list('thetas'=thetas))#,
                             # 'C_ZZe0'=C_ZZe0, 'C_ZY1'=C_ZY1,
                             # 'eigval_C_ZZe0'=eigval_C_ZZe0,
                             # 'eigval_C_eps'=eigval_C_eps,
                             # 'K'=K, 'K1'=K1,
                             # 'projector'=projector,
                             # 'projector_eps'=projector_eps,
                             # 'C_X'=C_X, 'C_eps'=C_eps) )
}


#' Moore–Penrose
#'
#' @param X
#' @param w
#' @param p
#'
#' @return
#' @export
#'
#' @examples
operatorARCHmp <- function(X, w=0.05, p=1){
  N <- ncol(X)
  n <- nrow(X)

  # Make tensor square data to use (n x n x N)
  X2_old <- array(apply(X, 2, tcrossprod),
                  dim=c(nrow(X),nrow(X),ncol(X)))
  X2 <- matrix(nrow=n^2,ncol=N) # (nn x N)
  for(i in 1:N){
    X2[,i] <- as.vector(X2_old[,,i])
  }

  # C_eps (n x n)
  C_eps <- matrix(NA, n, n)
  times <- 1:n/n
  for (i in 1:n){
    #C_eps[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
    C_eps[i,] <- exp(-abs(times[i]-times)/2)
  }
  # C_eps <- as.vector(C_eps) #(nn)

  #################
  ## Fit Model
  #################
  # C_X <- rowMeans(X2,dims=2)
  # # Define X2_eps
  # X2_eps <- X2
  # for(i in 1:N){
  #   X2_eps[,,i] <- X2[,,i] %*% C_eps / nrow(X2[,,i])
  # }
  #
  # # Define X2p
  # X2p <- X2p_eps <- array(dim=c(n*p,n,N))
  # for(k in 1:N){
  #   tmp <- tmp_eps <- list()
  #   for(j in 1:p){
  #     if(k-j+1 >0){ # Lag back if enough data
  #       tmp[[j]] <- X2[,,k-j+1]
  #       tmp_eps[[j]] <- X2_eps[,,k-j+1]
  #     } else{ # Can't lag back (too early)
  #       tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n, ncol=n)
  #     }
  #   }
  #   X2p[,,k] <- t(do.call("cbind", tmp))
  #   X2p_eps[,,k] <- t(do.call("cbind", tmp_eps))
  # }
  C_X <- rowMeans(X2)#,dims=2)
  # Define X2_eps
  X2_eps <- X2
  for(i in 1:N){
    # X2_eps[,,i] <- X2[,,i] %*% C_eps / nrow(X2[,,i])
    X2_eps[,i] <- as.vector(matrix(X2[,i], ncol=n,nrow=n) %*% C_eps) / n
  }

  # Define X2p
  #X2p <- X2p_eps <- array(dim=c(n*p,n,N))
  X2p <- X2p_eps <- matrix(nrow=p*n*n, ncol=N)
  for(k in 1:N){
    tmp <- tmp_eps <- list()
    for(j in 1:p){
      if(k-j+1 >0){ # Lag back if enough data
        # tmp[[j]] <- X2[,,k-j+1]
        # tmp_eps[[j]] <- X2_eps[,,k-j+1]
        tmp[[j]] <- X2[,k-j+1,drop=FALSE]
        tmp_eps[[j]] <- X2_eps[,k-j+1,drop=FALSE]
      } else{ # Can't lag back (too early)
        # tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n, ncol=n)
        tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n^2, ncol=1)
      }
    }
    # X2p[,,k] <- t(do.call("cbind", tmp))
    # X2p_eps[,,k] <- t(do.call("cbind", tmp_eps))
    X2p[,k] <- do.call("rbind", tmp)
    X2p_eps[,k] <- do.call("rbind", tmp_eps)
  }

  # mxp <- matrix(array(0,dim(X2p[,,1])))
  mxp <- matrix(0, ncol=1, nrow=nrow(X2p))
  for(l in p:N){
    # mxp <- mxp + X2p[,,l]
    mxp <- mxp + X2p[,l]
  }
  mxp <- mxp / (N-p+1)

  # Eqs 35 and 36
  tmp <- .estim_CD_hat(X = X2, Xp=X2p, Xpe = X2p_eps, p = p)
  C_hat <- tmp$C_hat
  D_hat <- tmp$D_hat

  # After Eq 36 and Eq 39
  eig_C <- eigen( C_hat %*% t(C_hat) )
  # eig_C <- eigen( C_hat )
  eig_Ce <- eigen( C_eps )
  CC_n <- nrow(eig_C$vectors)
  Ce_n <- nrow(eig_Ce$vectors)
  eig_C$vectors <- eig_C$vectors * sqrt(CC_n)
  eig_Ce$vectors <- eig_Ce$vectors * sqrt(Ce_n)

  K_C <- which.min(cumsum(Re(eig_C$values))/sum(Re(eig_C$values)) < 1-w)
  K_Ce <- which.min(cumsum(eig_Ce$values)/sum(eig_Ce$values) < 1-w)

  # Estimate alpha
  C_dag <- matrix(0, nrow=nrow(C_hat), ncol=ncol(C_hat))
  for(i in 1:K_C){
    C_dag <- C_dag +
      1 / eig_C$values[i] *
      ( eig_C$vectors[,i] %*% t(eig_C$vectors[,i]) )
  }
  alpha_hat <- apply(Re(D_hat %*% C_dag/K_C), 2, pmax, 0) # n (np) x np (np)

  ## Delta, Eq 39
  C_eps_dag <- matrix(0, nrow=nrow(C_eps), ncol=ncol(C_eps))
  for(i in 1:K_Ce){
    C_eps_dag <- C_eps_dag +
      1 / eig_Ce$values[i] *
      ( eig_Ce$vectors[,i] %*% t(eig_Ce$vectors[,i]) )
  }
  # alpha_sum <- matrix(0, nrow=nrow(C_X), ncol=ncol(C_X))
  # for(i in 1:p){
  #   alpha_sum <- alpha_sum +
  #     .apply_tensor(op=alpha_hat[,1:nrow(alpha_hat)+(i-1)*nrow(alpha_hat)],
  #                   X_apply=C_X)
  # }
  alpha_sum <- matrix(alpha_hat %*% mxp,
                      nrow = n,
                      ncol = n )
  delta_hat <-
    apply((C_X - (alpha_sum %*% C_eps) / n ) %*% C_eps_dag / K_Ce, 2, pmax, 0)

  list('hat_alpha'=alpha_hat,
       'hat_delta'=delta_hat,
       'error'=.compute_ARCH_error_op(
         X2 = X2, alpha = alpha_hat, delta = delta_hat,
         C_eps=C_eps, p=p, C_X=C_X))
}

###################################



.estimate_params_ARCH_op <- function(theta, vu, C_hat,D_hat, C_eps, C_X,
                                     projector_CC,projector_Ce, n, p, mxp){
  # Compute alpha, Eq 22
  C_hat_st <- t(Conj(C_hat)) # (pnn x pnn)
  C_dag <-
    C_hat_st %*% solve( C_hat %*% C_hat_st +
                          diag(theta,nrow(C_hat)) ) # (pnn x pnn)
  # (nn x pnn) x (pnn x pnn) x (pnn x pnn)
  alpha_hat <- (D_hat %*% C_dag %*% projector_CC)
  alpha_hat  <- apply(alpha_hat, 2, pmax, 0) # nn x nnp

  ## Delta, Eq 39
  C_eps_diag <- solve( C_eps + diag(vu, n) ) # n x n
  # alpha_sum <- .apply_tensor(alpha_hat, mxp) # (nn x nnp)x(np x n)
  alpha_sum <- matrix(alpha_hat %*% mxp,
                      nrow = sqrt(nrow(alpha_hat)),
                      ncol = sqrt(nrow(alpha_hat)) )

  delta_hat <-
    ((C_X - (alpha_sum %*% C_eps) / n ) %*%
       C_eps_diag  %*%
       projector_Ce)
  delta_hat <- apply(delta_hat, 2, pmax, 0)

  list('alpha_hat'=alpha_hat,
       'delta_hat'=delta_hat)
}

.estim_CD_hat <- function(X, Xp, Xpe, p){
  if(prod(dim(Xp) != dim(Xpe)))
    stop('Error: Xp and Xpe have different dimensions')

  # N <- dim(Xp)[3]
  N <- ncol(Xp)

  # Eqs 20 and 21
  #   C: (pnn x pnn),
  #   D: (nn x pnn),
  # Chat <-  matrix(0, nrow=nrow(Xp[,,1])*ncol(Xpe[,,1]),
  #                 ncol=nrow(Xp[,,1])*ncol(Xpe[,,1]))
  # Dhat <-  matrix(0, nrow=ncol(Xp[,,1])*nrow(X[,,1]),
  #                 ncol=nrow(Xp[,,1])*ncol(X[,,1]))
  Chat <-  matrix(0, nrow=nrow(Xp), ncol=nrow(Xp))
  Dhat <-  matrix(0, nrow=nrow(X), ncol=nrow(Xp))

  # After Eq 35
  mxp <- .mean_m(Xp,p=p,lag=0,lag2=0) # (100x20)
  mxpe <- .mean_m(Xpe,p=p,lag=0,lag2=0) # (100x20)
  mxppr <- .mean_m(Xp,p=p,lag=1,lag2=0) # (100x20)
  mxpr <- .mean_m(X,p=p,lag=1,lag2=1) # (20x20)

  # Chat <- tensorEVD::Kronecker(Xpe[,,N]-mxpe, t(Xp[,,N]-mxp))
  Chat <- (Xpe[,N]-mxpe) %*% t(Xp[,N]-mxp)
  for(k in p:(N-1)){
    ## C_hat, Eq 20
    #   N case computed before loop
    Chat <- Chat +
      # tensorEVD::Kronecker(Xpe[,,k]-mxpe,t(Xp[,,k]-mxp))
      (Xpe[,k]-mxpe) %*% t(Xp[,k]-mxp)

    ## D_hat, Eq 21
    Dhat <- Dhat +
      # tensorEVD::Kronecker(X[,,k+1]-mxpr, t(Xp[,,k]-mxppr))
      (X[,k+1]-mxpr) %*% t(Xp[,k]-mxppr)
  }

  list('C_hat'= Chat / (N-p),
       'D_hat'= Dhat / (N-p-1)
       )
}

.apply_tensor <- function(op,X_apply){
  matrix(op %*% matrixcalc::vec(X_apply)/ncol(op),
         nrow=nrow(X_apply),
         ncol=nrow(X_apply))
}

.mean_m <- function(X, p, lag,lag2=0){
  # dims <- dim(X)
  # mat <- matrix(0,nrow=dims[1],ncol=dims[2])
  mat <- matrix(0, nrow=nrow(X), ncol=1)
  N <- ncol(X)

  # for(k in p:(dims[3]-lag)){
  #   mat <- mat + X[,,k+lag2]
  # }
  for(k in p:(N-lag)){
    mat <- mat + X[,k+lag2]
  }

  mat / (N-p+1-lag)
}

.computeProjectorsAndInitialization <- function(C_hat, C_eps,w,n,p){
  # Between Eqs 22-23
  eig_CC <- eigen( C_hat %*% t(C_hat) ) # C for p=1
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
      ( eig_CC$vectors[,i] %*% t(eig_CC$vectors[,i]) )
  }
  for(i in 1:K_Ce){
    # (n x n)
    projector_Ce <- projector_Ce +
      ( eig_Ce$vectors[,i] %*% t(eig_Ce$vectors[,i]) )
  }

  ## Initialize
  theta <- eig_CC$value[1] / 2
  vu <- eig_Ce$value[1] / 2

  list(projector_CC=projector_CC,
       projector_Ce=projector_Ce,
       theta=theta,
       vu=vu)
}

##############################################

.compute_ARCH_error_op <- function(X2, alpha, delta, C_eps, C_X, p){
  # error <- rep(0, dim(X2)[3])
  error <- rep(0, ncol(X2))
  estimate <- .predict_ARCH_op(X2, alpha, delta, C_eps, p)

  # for(i in (1+p):dim(X2)[3]){
  for(i in (1+p):ncol(X2) ){
    # error[i] <- dot_integrate(apply((X2[,,i] - estimate[,,i])^2,
    #                               MARGIN = 2,
    #                               FUN = dot_integrate))
    error[i] <- dot_integrate((X2[,i] - estimate[,i])^2)
  }

  sum(error)
}

.predict_ARCH_op <- function(X2, alpha, delta, C_eps, p){
  # dims <- dim(X2)
  # Sigma <- array(NA, dims)
  N <- ncol(X2)
  n <- nrow(delta)
  Sigma <- array(NA, c(n,n,N))

  # Eq 2
  # for(i in (1+p):dims[3]){
  for(i in (1+p):N){
    Sigma[,,i] <- delta
    for(j in 1:p){
      Sigma[,,i] <- Sigma[,,i] +
        # .apply_tensor(op = alpha[,1:nrow(alpha)+(j-1)*nrow(alpha)],
        #               X_apply = X2[,,i-j])
        matrix(alpha[,1:nrow(alpha)+(j-1)*nrow(alpha)] %*% X2[,i-j],
               ncol=n,nrow=n)
    }
  }

  # Eq 3 (Eq 29)
  # Estimate <- array(NA, dims)
  Estimate <- matrix(NA, nrow=n^2, ncol=N)
  # for(i in (1+p):dims[3]){
  for(i in (1+p):N){
    Sigma_svd <- tryCatch({svd(Sigma[,,i])},
                          error=function(e){
                            svd::propack.svd(Sigma[,,i])
                          })
    Sigma_sqrt <- sqrt(n) * Sigma_svd$u %*%
      diag(sqrt(Sigma_svd$d)) %*%
      t(Sigma_svd$v)

    ## Eq 3
    # Estimate[,,i] <- (Sigma_sqrt %*% C_eps %*% Sigma_sqrt) / nrow(Sigma_sqrt)^2 #+ tmp
    Estimate[,i] <- (Sigma_sqrt %*% C_eps %*% Sigma_sqrt) / n^2
  }

  Estimate
}

