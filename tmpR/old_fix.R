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
operatorARCH4 <- function(X, w=0.05, thetas=NA, training=0.8, p=1){
  # Make tensor square data to use
  X2 <- array(apply(X, 2, tcrossprod),
              dim=c(nrow(X),nrow(X),ncol(X)))
  n <- nrow(X)

  # C_eps
  C_eps <- matrix(NA, n, n)
  times <- 1:n/n
  for (i in 1:n){
    C_eps[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
  }

  #################
  ## USE CV FOR MODEL SETUP
  #     Fit on part, test on pred
  #################
  if(training < 1){
    X_part_idx <- 1:(training * ncol(X))
    # lag back p so can test on the first new data!
    X_pred_idx <- (training * ncol(X)+1-p):ncol(X)
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
  X2p_part <- X2p_part_eps <- array(dim=c(n*p,n,N))
  for(k in 1:N){
    tmp <- tmp_eps <- list()
    for(j in 1:p){
      if(k-j+1 >0){ # We can look back at the lag
        tmp[[j]] <- X2_part[,,k-j+1]
        tmp_eps[[j]] <- X2_part_eps[,,k-j+1]
      } else{ # Too early in the sample to look back to data
        tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n, ncol=n)
      }
    }
    X2p_part[,,k] <- t(do.call("cbind", tmp))
    X2p_part_eps[,,k] <- t(do.call("cbind", tmp_eps))
  }

  # Eqs 35 and 36
  tmp <- .estim_CD_hat(X = X2_part, Xp=X2p_part, Xpe = X2p_part_eps, p = p)

  C_hat <- tmp$C_hat
  D_hat <- tmp$D_hat

  # After Eq 36 and Eq 39
  tmp <- .computeProjectors4(C_hat, C_eps, w, n,p)
  projector_CC <- tmp[[1]] # rows of C_hat
  projector_Ce <- tmp[[2]] # res of X2

  # Find the optimal Tikhonov parameters
  if(length(thetas)!=2 ||
     sum(is.na(thetas)) > 0 ){

    guess_theta <-
      eigen((C_hat %*% t(Conj(C_hat))) / nrow(C_hat))$value[1] / 2
    guess_vu <- eigen(C_eps)$value[1] / 2

    optim_results <-
      optim(c(guess_theta,guess_vu), fn=function(theta0, C_hat, D_hat,
                                                 projector_CC, projector_Ce,
                                                 n, p1, X21, C_X1){
        theta <- theta0[1]
        vu <- theta0[2]

        ests <- .estimate_params_ARCH_op3(theta = theta, vu=vu, C_hat = C_hat,
                                         D_hat = D_hat, C_eps=C_eps, C_X=C_X1,
                                         projector_CC = projector_CC,
                                         projector_Ce = projector_Ce,
                                         n = n, p = p1)

        .compute_ARCH_error_op(
          X2 = X21, alpha = ests$alpha_hat, delta = ests$delta_hat,
          C_eps=C_eps, p=p1, C_X=C_X1)
      },
      method="L-BFGS-B",#"Brent",
      lower=c(10^-8,0),upper=c(Inf,Inf),
      C_hat=C_hat, D_hat=D_hat, projector_CC=projector_CC,
      projector_Ce=projector_Ce, n=n, p1=p, X21=X2_pred, C_X1=C_part_X)

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

  # Define X2p
  X2p <- X2p_eps <- array(dim=c(n*p,n,N))
  for(k in 1:N){
    tmp <- tmp_eps <- list()
    for(j in 1:p){
      if(k-j+1 >0){ # Lag back if enough data
        tmp[[j]] <- X2[,,k-j+1]
        tmp_eps[[j]] <- X2_eps[,,k-j+1]
      } else{ # Can't lag back (too early)
        tmp[[j]] <- tmp_eps[[j]] <- matrix(0, nrow=n, ncol=n)
      }
    }
    X2p[,,k] <- t(do.call("cbind", tmp))
    X2p_eps[,,k] <- t(do.call("cbind", tmp_eps))
  }

  # Eqs 35 and 36
  tmp <- .estim_CD_hat(X = X2, Xp=X2p, Xpe = X2p_eps, p = p)

  C_hat <- tmp$C_hat
  D_hat <- tmp$D_hat

  tmp <- .computeProjectors4(C_hat, C_eps,w,n,p)
  projector_CC <- tmp[[1]]
  projector_Ce <- tmp[[2]]

  ests <- .estimate_params_ARCH_op3(theta = thetas[1], vu=thetas[2],
                                   C_hat = C_hat,
                                   D_hat = D_hat, C_eps=C_eps, C_X=C_X,
                                   projector_CC = projector_CC,
                                   projector_Ce = projector_Ce,
                                   n = n, p = p)

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

.computeProjectors4 <- function(C_hat, C_eps,w,n,p){
  # After Eqs 36 and 39
  eig_CC <- eigen( C_hat %*% t(C_hat) )
  eig_Ce <- eigen( C_eps )
  CC_n <- nrow(eig_CC$vectors)
  Ce_n <- nrow(eig_Ce$vectors)
  eig_CC$vectors <- eig_CC$vectors * sqrt(CC_n)
  eig_Ce$vectors <- eig_Ce$vectors * sqrt(Ce_n)

  K_CC <- which.min(cumsum(eig_CC$values)/sum(eig_CC$values) < 1-w)
  K_Ce <- which.min(cumsum(eig_Ce$values)/sum(eig_Ce$values) < 1-w)

  # Compute Projection
  projector_CC <- matrix(0, nrow=CC_n, ncol=CC_n)
  projector_Ce <- matrix(0, nrow=Ce_n, ncol=Ce_n)

  for(i in 1:K_CC){
    projector_CC <- projector_CC +
      ( eig_CC$vectors[,i] %*% t(eig_CC$vectors[,i]) )
  }
  for(i in 1:K_Ce){
    projector_Ce <- projector_Ce +
      ( eig_Ce$vectors[,i] %*% t(eig_Ce$vectors[,i]) )
  }

  list(projector_CC, projector_Ce)
}


.estimate_params_ARCH_op3 <- function(theta, vu, C_hat,D_hat, C_eps, C_X,
                                      projector_CC,projector_Ce, n, p){
  # Compute alpha, Eq 36
  C_hat_st <- t(Conj(C_hat))
  C_dag <-
    C_hat_st %*% solve( C_hat %*% C_hat_st +
                          diag(theta,nrow(C_hat)) )
  alpha_hat <- (D_hat %*% C_dag %*% projector_CC)
  alpha_hat  <- apply(alpha_hat, 2, pmax, 0) # n (np) x np (np)

  ## Delta, Eq 39
  C_eps_diag <- solve( C_eps + diag(vu, n) )
  alpha_sum <- matrix(0, nrow=nrow(C_X), ncol=ncol(C_X))
  for(i in 1:p){
    alpha_sum <- alpha_sum +
      .apply_tensor(op=alpha_hat[,1:nrow(alpha_hat)+(i-1)*nrow(alpha_hat)],
                    X_apply=C_X)
  }
  delta_hat <-
    ((C_X - (alpha_sum %*% C_eps) / n ) %*%
       C_eps_diag  %*%
       projector_Ce)
  delta_hat <- apply(delta_hat, 2, pmax, 0)

  list('alpha_hat'=alpha_hat,
       'delta_hat'=delta_hat)
}
