#' # data <- generateFunctionalCondVol(21,500,'arch')
#' # res <- operatorARCH(X=data$garch_mat,p=1)
#' # or <- .get_oracle(data)
#' #
#' # data <- generateFunctionalCondVol_ff(21,1000,'arch')#,delta_par = function(s,t){min(s,t)})
#' # res <- operatorARCH(X=data$garch_mat,p=2)
#' # res <- operatorARCH(X=data$garch_mat,p=1)
#' # or <- .get_oracle_op(data)
#' #
#' # plot_fd(res$hat_alpha)
#' # plot_fd(or$hat_alpha)
#' #
#' # plot_fd(res$hat_delta)
#' # plot_fd(or$hat_delta)
NULL

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
operatorGARCH <- function(X, w=0.05, thetas=NA, training=0.8, L=10){
  # Make tensor square data to use
  X2 <- array(apply(X, 2, tcrossprod),
              dim=c(nrow(X),nrow(X),ncol(X)))
  n <- nrow(X)

  # C_eps
  C_eps <- matrix(NA, n, n)
  times <- 1:n/n
  for (i in 1:n){
    #C_eps[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
    C_eps[i,] <- exp(-abs(times[i]-times)/2)
  }

  #################
  ## USE CV FOR MODEL SETUP
  #     Fit on part, test on pred
  #################
  if(training < 1){
    X_part_idx <- 1:(training * ncol(X))
    # lag back p so can test on the first new data!
    X_pred_idx <- (training * ncol(X)+1-1):ncol(X)
  }else{
    X_part_idx <- 1:ncol(X)
    X_pred_idx <- 1:ncol(X)
  }
  N <- length(X_part_idx)
  C_part_X <- rowMeans(X2[,,X_part_idx],dims=2)

  X2_part <- X2[,,X_part_idx]
  X2_pred <- X2[,,X_pred_idx]

  # Define X2_eps (after Eq 29)
  X2_part_eps <- X2_part
  for(i in 1:N){
    X2_part_eps[,,i] <- (X2_part[,,i] %*% C_eps) / n
  }

  # Define X2p and X2p_eps (Eqs 30 and 31)
  X2L_part <- X2L_part_eps <- array(dim=c(n*L,n,N))
  for(k in 1:N){
    tmp <- tmp_eps <- list()
    for(j in 1:L){
      if(k-j+1 >0){ # We can look back at the lag
        tmp[[j]] <- X2_part[,,k-j+1]
        tmp_eps[[j]] <- X2_part_eps[,,k-j+1]
      } else{ # Too early in the sample to look back to data
        tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n, ncol=n)
      }
    }
    X2L_part[,,k] <- t(do.call("cbind", tmp))
    X2L_part_eps[,,k] <- t(do.call("cbind", tmp_eps))
  }

  # Eqs 35 and 36
  tmp <- .estim_CD_hat(X = X2_part, Xp=X2L_part, Xpe = X2L_part_eps, p = L)

  C_hat <- tmp$C_hat
  D_hat <- tmp$D_hat

  # After Eq 36 and Eq 39
  tmp <- .computeProjectorsAndInitialization(C_hat, C_eps, w, n, L)
  projector_CC <- tmp$projector_CC # rows of C_hat
  projector_Ce <- tmp$projector_Ce # res of X2
  guess_theta <- tmp$theta
  guess_vu <- tmp$vu

  # Find the optimal Tikhonov parameters
  if(length(thetas)!=3 ||
     sum(is.na(thetas)) > 0 ){

    guess_gamma <- .getGammaEst(C_hat,D_hat,projector_CC,guess_theta,n)

    # library("optimParallel")
    # cl <- makeCluster(4)     # set the number of processor cores
    # setDefaultCluster(cl=cl) # set 'cl' as default cluster
    # optimParallel(c(guess_theta,guess_gamma,guess_vu), #fn=negll, x=x,
                  #method = "L-BFGS-B", lower=c(-Inf, .0001))
    optim_results <-
    optim(c(guess_theta,guess_gamma,guess_vu),
            fn=function(theta0, C_hat, D_hat,
                        projector_CC, projector_Ce,
                        n, L1, X21, C_X1,w1){
              # devtools::load_all()
              theta <- theta0[1]
              gamma <- theta0[2]
              vu <- theta0[3]

              ests <- .estimate_params_GARCH_op(
                theta = theta, gamma=gamma, vu=vu,
                C_hat = C_hat, D_hat = D_hat, C_eps=C_eps, C_X=C_X1,
                projector_CC = projector_CC, projector_Ce = projector_Ce,
                n = n, L = L1, w=w1)

              .compute_GARCH_error_op(
                X2 = X21, alpha = ests$alpha_hat,
                beta = ests$beta_hat, delta = ests$delta_hat,
                C_eps=C_eps)
            },
            method="L-BFGS-B",#"Brent",
            lower=c(10^-8,10^-8,0),upper=c(Inf,Inf,Inf),
            C_hat=C_hat, D_hat=D_hat, projector_CC=projector_CC,
            projector_Ce=projector_Ce, n=n, L1=L, X21=X2_pred,
            C_X1=C_part_X, w1=w)

    thetas <- optim_results$par
  }

  #################
  ## USE ALL DATA FOR FINAL SETUP
  #################

  N <- ncol(X) # Number of observations
  C_X <- rowMeans(X2,dims=2)
  # Define X2_eps
  X2_eps <- X2
  for(i in 1:N){
    X2_eps[,,i] <- X2[,,i] %*% C_eps / nrow(X2[,,i])
  }

  # Define X2L
  X2L <- X2L_eps <- array(dim=c(n*L,n,N))
  for(k in 1:N){
    tmp <- tmp_eps <- list()
    for(j in 1:L){
      if(k-j+1 >0){ # Lag back if enough data
        tmp[[j]] <- X2[,,k-j+1]
        tmp_eps[[j]] <- X2_eps[,,k-j+1]
      } else{ # Can't lag back (too early)
        tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n, ncol=n)
      }
    }
    X2L[,,k] <- t(do.call("cbind", tmp))
    X2L_eps[,,k] <- t(do.call("cbind", tmp_eps))
  }

  # Eqs 35 and 36
  tmp <- .estim_CD_hat(X = X2, Xp=X2L, Xpe = X2L_eps, p = L)

  C_hat <- tmp$C_hat
  D_hat <- tmp$D_hat

  tmp <- .computeProjectorsAndInitialization(C_hat, C_eps, w, n, L)
  projector_CC <- tmp$projector_CC # rows of C_hat
  projector_Ce <- tmp$projector_Ce # res of X2

  ests <- .estimate_params_GARCH_op(
    theta = thetas[1],gamma=thetas[2], vu=thetas[3],
    C_hat = C_hat, D_hat = D_hat, C_eps=C_eps, C_X=C_X,
    projector_CC = projector_CC, projector_Ce = projector_Ce,
    n = n, L = L, w=w)

  list('hat_alpha'=ests$alpha_hat,
       'hat_beta'=ests$beta_hat,
       'hat_delta'=ests$delta_hat,
       'error'=.compute_GARCH_error_op(
         X2 = X2, alpha = ests$alpha_hat,
         beta = ests$beta_hat, delta = ests$delta_hat,
         C_eps=C_eps),
       'fitting_params'=list('thetas'=thetas) )
}

operatorGARCH1 <- function(X, w=0.05, thetas=NA, training=0.8, L=10){
  # Make tensor square data to use
  X2 <- array(apply(X, 2, tcrossprod),
              dim=c(nrow(X),nrow(X),ncol(X)))
  n <- nrow(X)

  # C_eps
  C_eps <- matrix(NA, n, n)
  times <- 1:n/n
  for (i in 1:n){
    #C_eps[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
    C_eps[i,] <- exp(-abs(times[i]-times)/2)
  }

  #################
  ## USE CV FOR MODEL SETUP
  #     Fit on part, test on pred
  #################
  if(training < 1){
    X_part_idx <- 1:(training * ncol(X))
    # lag back p so can test on the first new data!
    X_pred_idx <- (training * ncol(X)+1-1):ncol(X)
  }else{
    X_part_idx <- 1:ncol(X)
    X_pred_idx <- 1:ncol(X)
  }
  N <- length(X_part_idx)
  C_part_X <- rowMeans(X2[,,X_part_idx],dims=2)

  X2_part <- X2[,,X_part_idx]
  X2_pred <- X2[,,X_pred_idx]

  # Define X2_eps (after Eq 29)
  X2_part_eps <- X2_part
  for(i in 1:N){
    X2_part_eps[,,i] <- (X2_part[,,i] %*% C_eps) / n
  }

  # Define X2p and X2p_eps (Eqs 30 and 31)
  X2L_part <- X2L_part_eps <- array(dim=c(n*L,n,N))
  for(k in 1:N){
    tmp <- tmp_eps <- list()
    for(j in 1:L){
      if(k-j+1 >0){ # We can look back at the lag
        tmp[[j]] <- X2_part[,,k-j+1]
        tmp_eps[[j]] <- X2_part_eps[,,k-j+1]
      } else{ # Too early in the sample to look back to data
        tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n, ncol=n)
      }
    }
    X2L_part[,,k] <- t(do.call("cbind", tmp))
    X2L_part_eps[,,k] <- t(do.call("cbind", tmp_eps))
  }

  # Eqs 35 and 36
  tmp <- .estim_CD_hat(X = X2_part, Xp=X2L_part, Xpe = X2L_part_eps, p = L)

  C_hat <- tmp$C_hat
  D_hat <- tmp$D_hat

  # After Eq 36 and Eq 39
  tmp <- .computeProjectors(C_hat, C_eps, w, n, L)
  projector_CC <- tmp[[1]] # rows of C_hat
  projector_Ce <- tmp[[2]] # res of X2

  # Find the optimal Tikhonov parameters
  if(length(thetas)!=3 ||
     sum(is.na(thetas)) > 0 ){

    guess_theta <-
      eigen((C_hat %*% t(Conj(C_hat))) / nrow(C_hat))$value[1] / 2
    guess_gamma <- .getGammaEst(C_hat,D_hat,projector_CC,guess_theta,n)
    guess_vu <- eigen(C_eps)$value[1] / 2

    optim_results <-
      optim(c(guess_theta,guess_gamma,guess_vu),
            fn=function(theta0, C_hat, D_hat,
                        projector_CC, projector_Ce,
                        n, L1, X21, C_X1,w1){
              theta <- theta0[1]
              gamma <- theta0[2]
              vu <- theta0[3]

              ests <- .estimate_params_GARCH_op1(
                theta = theta, gamma=gamma, vu=vu,
                C_hat = C_hat, D_hat = D_hat, C_eps=C_eps, C_X=C_X1,
                projector_CC = projector_CC, projector_Ce = projector_Ce,
                n = n, L = L1, w=w1)

              .compute_GARCH_error_op(
                X2 = X21, alpha = ests$alpha_hat,
                beta = ests$beta_hat, delta = ests$delta_hat,
                C_eps=C_eps)
            },
            method="L-BFGS-B",#"Brent",
            lower=c(10^-8,10^-8,0),upper=c(Inf,Inf,Inf),
            C_hat=C_hat, D_hat=D_hat, projector_CC=projector_CC,
            projector_Ce=projector_Ce, n=n, L1=L, X21=X2_pred, C_X1=C_part_X, w1=w)

    thetas <- optim_results$par
  }

  #################
  ## USE ALL DATA FOR FINAL SETUP
  #################

  N <- ncol(X) # Number of observations
  C_X <- rowMeans(X2,dims=2)
  # Define X2_eps
  X2_eps <- X2
  for(i in 1:N){
    X2_eps[,,i] <- X2[,,i] %*% C_eps / nrow(X2[,,i])
  }

  # Define X2L
  X2L <- X2L_eps <- array(dim=c(n*L,n,N))
  for(k in 1:N){
    tmp <- tmp_eps <- list()
    for(j in 1:L){
      if(k-j+1 >0){ # Lag back if enough data
        tmp[[j]] <- X2[,,k-j+1]
        tmp_eps[[j]] <- X2_eps[,,k-j+1]
      } else{ # Can't lag back (too early)
        tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n, ncol=n)
      }
    }
    X2L[,,k] <- t(do.call("cbind", tmp))
    X2L_eps[,,k] <- t(do.call("cbind", tmp_eps))
  }

  # Eqs 35 and 36
  tmp <- .estim_CD_hat(X = X2, Xp=X2L, Xpe = X2L_eps, p = L)

  C_hat <- tmp$C_hat
  D_hat <- tmp$D_hat

  tmp <- .computeProjectors(C_hat, C_eps,w,n,L)
  projector_CC <- tmp[[1]]
  projector_Ce <- tmp[[2]]

  ests <- .estimate_params_GARCH_op1(
    theta = thetas[1],gamma=thetas[2], vu=thetas[3],
    C_hat = C_hat, D_hat = D_hat, C_eps=C_eps, C_X=C_X,
    projector_CC = projector_CC, projector_Ce = projector_Ce,
    n = n, L = L, w=w)

  list('hat_alpha'=ests$alpha_hat,
       'hat_beta'=ests$beta_hat,
       'hat_delta'=ests$delta_hat,
       'error'=.compute_GARCH_error_op(
         X2 = X2, alpha = ests$alpha_hat,
         beta = ests$beta_hat, delta = ests$delta_hat,
         C_eps=C_eps),
       'fitting_params'=list('thetas'=thetas) )
}

#' Title
#'
#' @param X
#' @param w
#' @param thetas
#' @param training
#' @param L
#'
#' @return
#' @export
#'
#' @examples
operatorGARCHmp <- function(X, w=0.05, thetas=NA, training=0.8, L=5){
  # Make tensor square data to use
  X2 <- array(apply(X, 2, tcrossprod),
              dim=c(nrow(X),nrow(X),ncol(X)))
  n <- nrow(X)

  # C_eps
  C_eps <- matrix(NA, n, n)
  times <- 1:n/n
  for (i in 1:n){
    #C_eps[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
    C_eps[i,] <- exp(-abs(times[i]-times)/2)
  }

  #################
  ## Fit Model
  #################
  N <- ncol(X) # Number of observations
  C_X <- rowMeans(X2,dims=2)
  # Define X2_eps
  X2_eps <- X2
  for(i in 1:N){
    X2_eps[,,i] <- X2[,,i] %*% C_eps / nrow(X2[,,i])
  }

  # Define X2L
  X2L <- X2L_eps <- array(dim=c(n*L,n,N))
  for(k in 1:N){
    tmp <- tmp_eps <- list()
    for(j in 1:L){
      if(k-j+1 >0){ # Lag back if enough data
        tmp[[j]] <- X2[,,k-j+1]
        tmp_eps[[j]] <- X2_eps[,,k-j+1]
      } else{ # Can't lag back (too early)
        tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n, ncol=n)
      }
    }
    X2L[,,k] <- t(do.call("cbind", tmp))
    X2L_eps[,,k] <- t(do.call("cbind", tmp_eps))
  }

  # Eqs 35 and 36
  tmp <- .estim_CD_hat(X = X2, Xp=X2L, Xpe = X2L_eps, p = L)

  C_hat <- tmp$C_hat
  D_hat <- tmp$D_hat

  # After Eq 36 and Eq 39
  eig_C <- eigen( C_hat )
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
  xi <- Re(D_hat %*% C_dag)
  xi1 <- xi[,1:nrow(xi)]
  xi2 <- xi[,nrow(xi) + 1:nrow(xi)]

  alpha_hat  <- apply(xi1, 2, pmax, 0) # n (np) x np (np),p=1

  ## Beta, Eq 26
  eig_ei1 <- eigen(xi1)
  ei1_n <- nrow(eig_ei1$vectors)
  eig_ei1$vectors <- eig_ei1$vectors * sqrt(ei1_n)
  K_ei1 <- which.min(cumsum(Re(eig_ei1$values))/sum(Re(eig_ei1$values)) < 1-w)

  C_dag1 <- matrix(0, nrow=nrow(eig_ei1$vectors), ncol=ncol(eig_ei1$vectors))
  for(i in 1:K_ei1){
    C_dag1 <- C_dag1 +
      1 / eig_ei1$values[i] *
      ( eig_ei1$vectors[,i] %*% t(eig_ei1$vectors[,i]) )
  }

  beta_hat <- apply( Re(xi2 %*% C_dag1), 2, pmax, 0)


  ## Delta, Eq 28
  C_eps_dag <- matrix(0, nrow=nrow(C_eps), ncol=ncol(C_eps))
  for(i in 1:K_C){
    C_eps_dag <- C_eps_dag +
      1 / eig_Ce$values[i] *
      ( eig_Ce$vectors[,i] %*% t(eig_Ce$vectors[,i]) )
  }


  alpha_sum <- .apply_tensor(op=alpha_hat, X_apply=C_X)
  beta_sum <-
    .apply_tensor(op=beta_hat,
                  X_apply=C_X %*% C_eps_dag)

  delta_hat <-
    (C_X - alpha_sum %*% C_eps - beta_sum %*% C_eps) %*%
    C_eps_dag
  delta_hat <- apply(delta_hat, 2, pmax, 0)

  list('hat_alpha'=alpha_hat,
       'hat_beta'=beta_hat,
       'hat_delta'=delta_hat,
       'error'=.compute_GARCH_error_op(
         X2 = X2, alpha = alpha_hat,
         beta = beta_hat, delta = delta_hat,
         C_eps=C_eps) )
}

operatorGARCHmp1 <- function(X, w=0.05, thetas=NA, training=0.8, L=10){
  # Make tensor square data to use
  X2 <- array(apply(X, 2, tcrossprod),
              dim=c(nrow(X),nrow(X),ncol(X)))
  n <- nrow(X)

  # C_eps
  C_eps <- matrix(NA, n, n)
  times <- 1:n/n
  for (i in 1:n){
    #C_eps[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
    C_eps[i,] <- exp(-abs(times[i]-times)/2)
  }

  #################
  ## Fit Model
  #################
  N <- ncol(X) # Number of observations
  C_X <- rowMeans(X2,dims=2)
  # Define X2_eps
  X2_eps <- X2
  for(i in 1:N){
    X2_eps[,,i] <- X2[,,i] %*% C_eps / nrow(X2[,,i])
  }

  # Define X2L
  X2L <- X2L_eps <- array(dim=c(n*L,n,N))
  for(k in 1:N){
    tmp <- tmp_eps <- list()
    for(j in 1:L){
      if(k-j+1 >0){ # Lag back if enough data
        tmp[[j]] <- X2[,,k-j+1]
        tmp_eps[[j]] <- X2_eps[,,k-j+1]
      } else{ # Can't lag back (too early)
        tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n, ncol=n)
      }
    }
    X2L[,,k] <- t(do.call("cbind", tmp))
    X2L_eps[,,k] <- t(do.call("cbind", tmp_eps))
  }

  # Eqs 35 and 36
  tmp <- .estim_CD_hat(X = X2, Xp=X2L, Xpe = X2L_eps, p = L)

  C_hat <- tmp$C_hat
  D_hat <- tmp$D_hat

  # After Eq 36 and Eq 39
  eig_C <- eigen( C_hat )
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
  xi <- Re(D_hat %*% C_dag)
  xi1 <- xi[,1:nrow(xi)]
  xi2 <- xi[,nrow(xi) + 1:nrow(xi)]

  alpha_hat  <- apply(xi1, 2, pmax, 0) # n (np) x np (np),p=1

  ## Beta, Eq 26
  eig_ei1 <- eigen(xi1)
  ei1_n <- nrow(eig_ei1$vectors)
  eig_ei1$vectors <- eig_ei1$vectors * sqrt(ei1_n)
  K_ei1 <- which.min(cumsum(Re(eig_ei1$values))/sum(Re(eig_ei1$values)) < 1-w)

  C_dag1 <- matrix(0, nrow=nrow(eig_ei1$vectors), ncol=ncol(eig_ei1$vectors))
  for(i in 1:K_ei1){
    C_dag1 <- C_dag1 +
      1 / eig_ei1$values[i] *
      ( eig_ei1$vectors[,i] %*% t(eig_ei1$vectors[,i]) )
  }

  beta_hat <- apply( Re(xi2 %*% C_dag1), 2, pmax, 0)


  ## Delta, Eq 28
  C_eps_dag <- matrix(0, nrow=nrow(C_eps), ncol=ncol(C_eps))
  for(i in 1:K_C){
    C_eps_dag <- C_eps_dag +
      1 / eig_Ce$values[i] *
      ( eig_Ce$vectors[,i] %*% t(eig_Ce$vectors[,i]) )
  }


  alpha_sum <- .apply_tensor(op=alpha_hat, X_apply=C_X)
  beta_sum <-
    .apply_tensor(op=beta_hat,
                  X_apply=C_X %*% C_eps_dag)

  delta_hat <-
    (C_X - alpha_sum %*% C_eps/n - beta_sum %*% C_eps/n) %*%
    C_eps_dag
  delta_hat <- apply(delta_hat, 2, pmax, 0)

  list('hat_alpha'=alpha_hat,
       'hat_beta'=beta_hat,
       'hat_delta'=delta_hat,
       'error'=.compute_GARCH_error_op(
         X2 = X2, alpha = alpha_hat,
         beta = beta_hat, delta = delta_hat,
         C_eps=C_eps) )
}

######################################

.getGammaEst <- function(C_hat,D_hat,projector_CC,theta,n){

  C_hat_st <- t(Conj(C_hat))
  C_dag <- C_hat_st %*% solve( C_hat %*% C_hat_st + diag(theta,nrow(C_hat)) )

  xi <- (D_hat %*% C_dag %*% projector_CC)
  xi1 <- xi[,1:nrow(xi)]

  eigen( xi1 %*% t(xi1) )$values[1] / 2
}

.estimate_params_GARCH_op <- function(theta, gamma, vu, C_hat,D_hat, C_eps, C_X,
                                      projector_CC,projector_Ce, n, L, w){
  # Compute alpha, Eq 26
  C_hat_st <- t(Conj(C_hat))
  C_dag <- C_hat_st %*% solve( C_hat %*% C_hat_st + diag(theta,nrow(C_hat)) )

  xi <- (D_hat %*% C_dag %*% projector_CC)
  xi1 <- xi[,1:nrow(xi)]
  xi2 <- xi[,nrow(xi) + 1:nrow(xi)]

  alpha_hat <- xi1
  alpha_hat  <- apply(alpha_hat, 2, pmax, 0) # n (np) x np (np),p=1

  ## Beta, Eq 26
  xi1_st <- t(Conj(xi1))
  xi1_dag <- xi1_st %*% solve( xi1 %*% xi1_st + diag(gamma,nrow(xi1)) )
  projector_xi1xi1 <- .computeProjectorxi1(xi1,w)

  beta_hat <- xi2 %*% xi1_dag %*% projector_xi1xi1
  beta_hat <- apply(beta_hat, 2, pmax, 0)

  ## Delta, Eq 28
  C_eps_diag <- solve( C_eps + diag(vu, n) )
  #alpha_sum <- matrix(0, nrow=nrow(C_X), ncol=ncol(C_X))
  alpha_sum <- .apply_tensor(op=alpha_hat, X_apply=C_X)
  beta_sum <-
    .apply_tensor(op=beta_hat,
                  X_apply=C_X %*% C_eps_diag %*% projector_Ce) #TODO:: CHECK

  delta_hat <-
    (C_X - (alpha_sum - beta_sum) %*% C_eps / n) %*% # TODO:: CHeck
    C_eps_diag  %*%
    projector_Ce
  delta_hat <- apply(delta_hat, 2, pmax, 0)

  list('alpha_hat'=alpha_hat,
       'beta_hat'=beta_hat,
       'delta_hat'=delta_hat)
}

.estimate_params_GARCH_op1 <- function(theta, gamma, vu, C_hat,D_hat, C_eps, C_X,
                                      projector_CC,projector_Ce, n, L, w){
  # Compute alpha, Eq 26
  C_hat_st <- t(Conj(C_hat))
  C_dag <- C_hat_st %*% solve( C_hat %*% C_hat_st + diag(theta,nrow(C_hat)) )

  xi <- (D_hat %*% C_dag %*% projector_CC)
  xi1 <- xi[,1:nrow(xi)]
  xi2 <- xi[,nrow(xi) + 1:nrow(xi)]

  alpha_hat <- xi1
  alpha_hat  <- apply(alpha_hat, 2, pmax, 0) # n (np) x np (np),p=1

  ## Beta, Eq 26
  xi1_st <- t(Conj(xi1))
  xi1_dag <- xi1_st %*% solve( xi1 %*% xi1_st + diag(gamma,nrow(xi1)) )
  projector_xi1xi1 <- .computeProjectorxi1(xi1,w)

  beta_hat <- xi2 %*% xi1_dag %*% projector_xi1xi1
  beta_hat <- apply(beta_hat, 2, pmax, 0)

  ## Delta, Eq 28
  C_eps_diag <- solve( C_eps + diag(vu, n) )
  #alpha_sum <- matrix(0, nrow=nrow(C_X), ncol=ncol(C_X))
  alpha_sum <- .apply_tensor(op=alpha_hat, X_apply=C_X)
  beta_sum <-
    .apply_tensor(op=beta_hat,
                  X_apply=C_X %*% C_eps_diag %*% projector_Ce)

  delta_hat <-
    (C_X - alpha_sum %*% C_eps/n - beta_sum %*% C_eps/n) %*%
    C_eps_diag  %*%
    projector_Ce
  delta_hat <- apply(delta_hat, 2, pmax, 0)

  list('alpha_hat'=alpha_hat,
       'beta_hat'=beta_hat,
       'delta_hat'=delta_hat)
}

.computeProjectorxi1 <- function(xi1,w){
  eig_xi1 <- eigen( xi1 %*% t(xi1) )
  xi1_n <- nrow(eig_xi1$vectors)
  eig_xi1$vectors <- eig_xi1$vectors * sqrt(xi1_n)
  K_xi1 <- which.min(cumsum(eig_xi1$values)/sum(eig_xi1$values) < 1-w)

  # Compute Projection
  projector_xi1xi1 <- matrix(0, nrow=xi1_n, ncol=xi1_n)

  for(i in 1:K_xi1){
    projector_xi1xi1 <- projector_xi1xi1 +
      ( eig_xi1$vectors[,i] %*% t(eig_xi1$vectors[,i]) )
  }

  projector_xi1xi1
}

##############################################

.compute_GARCH_error_op <- function(X2, alpha, beta, delta, C_eps){
  error <- rep(0, dim(X2)[3])
  estimate <- .predict_GARCH_op(X2, alpha, beta, delta, C_eps)

  for(i in 2:dim(X2)[3]){
    error[i] <- dot_integrate(apply((X2[,,i] - estimate[,,i])^2,
                                    MARGIN = 2,
                                    FUN = dot_integrate))
  }

  sum(error)
}

.predict_GARCH_op <- function(X2, alpha, beta, delta, C_eps){
  dims <- dim(X2)
  Sigma <- array(NA, dims)

  # Eq 2
  Sigma[,,1] <-  array(0, dims[1:2])
  for(i in 2:dims[3]){
    Sigma[,,i] <- delta +
      .apply_tensor(op = alpha,
                    X_apply = X2[,,i-1]) +
      .apply_tensor(op = beta,
                    X_apply = Sigma[,,i-1])
                    #X_apply = X2[,,i-1] %*% solve(C_eps))
  }

  # Eq 3 (Eq 29)
  Estimate <- array(NA, dims)
  for(i in 2:dims[3]){
    Sigma_svd <- tryCatch({svd(Sigma[,,i])},
                          error=function(e){
                            svd::propack.svd(Sigma[,,i])
                          })
    Sigma_sqrt <- sqrt(nrow(X2)) * Sigma_svd$u %*%
      diag(sqrt(Sigma_svd$d)) %*%
      t(Sigma_svd$v)

    ## Eq 3
    Estimate[,,i] <- (Sigma_sqrt %*% C_eps %*% Sigma_sqrt) / nrow(Sigma_sqrt)^2
  }

  Estimate
}

