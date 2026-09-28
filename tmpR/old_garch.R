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
  if(length(thetas)!=2 ||
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

        ests <- .estimate_params_GARCH_op(
          theta = theta, gamma=gamma, vu=vu,
          C_hat = C_hat, _hat = D_hat, C_eps=C_eps, C_X=C_X1,
          projector_CC = projector_CC, projector_Ce = projector_Ce,
          n = n, L = L1, w=w1)

        .compute_GARCH_error_op(
          X2 = X21, alpha = ests$alpha_hat,
          beta = ests$beta_hat, delta = ests$delta_hat,
          C_eps=C_eps, L=L1, C_X=C_X1)
      },
      method="L-BFGS-B",#"Brent",
      lower=c(10^-8,0),upper=c(Inf,Inf),
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
  tmp <- .estim_CD_hat(X = X2, XL=X2L, XLe = X2L_eps, L = L)

  C_hat <- tmp$C_hat
  D_hat <- tmp$D_hat

  tmp <- .computeProjectors(C_hat, C_eps,w,n,L)
  projector_CC <- tmp[[1]]
  projector_Ce <- tmp[[2]]

  ests <- .estimate_params_GARCH_op(
    theta = thetas[1],gamma=thetas[2], vu=thetas[3],
    C_hat = C_hat, D_hat = D_hat, C_eps=C_eps, C_X=C_X,
    projector_CC = projector_CC, projector_Ce = projector_Ce,
    n = n, L = L, w=w)

  list('hat_alpha'=ests$alpha_hat,
       'hat_beta'=ests$beta_hat,
       'hat_delta'=ests$delta_hat,
       'error'=.compute_GARCH_error_op(
         X2 = X21, alpha = ests$alpha_hat,
         beta = ests$beta_hat, delta = ests$delta_hat,
         C_eps=C_eps, L=L1, C_X=C_X1),
       'fitting_params'=list('thetas'=thetas) )
}

.getGammaEst <- function(C_hat,D_hat,projector_CC,theta,n){

  C_hat_st <- t(Conj(C_hat))
  C_dag <- C_hat_st %*% solve( C_hat %*% C_hat_st + diag(theta,n) )

  xi <- (D_hat %*% C_dag %*% projector_CC)
  xi1 <- xi[,1:nrow(xi)]

  eigen( xi1 %*% t(xi1) )[1] / 2
}

.estimate_params_GARCH_op <- function(theta, gamma, vu, C_hat,D_hat, C_eps, C_X,
                                     projector_CC,projector_Ce, n, L, w){
  # Compute alpha, Eq 26
  C_hat_st <- t(Conj(C_hat))
  C_dag <- C_hat_st %*% solve( C_hat %*% C_hat_st + diag(theta,n) )

  xi <- (D_hat %*% C_dag %*% projector_CC)
  xi1 <- xi[,1:nrow(xi)]
  xi2 <- xi[,nrow(xi) + 1:nrow(xi)]

  alpha_hat <- xi1
  alpha_hat  <- apply(alpha_hat, 2, pmax, 0) # n (np) x np (np),p=1

  ## Beta, Eq 26
  xi1_st <- t(Conj(xi1))
  xi1_dag <- xi1_st %*% solve( xi1 %*% xi1_st + diag(gamma,n) )
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
    (C_X - alpha_sum %*% C_eps - beta_hat %*% C_eps) %*%
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
  projector_xi1xi1 <- matrix(0, nrow=CC_n, ncol=CC_n)

  for(i in 1:K_xi1){
    projector_xi1xi1 <- projector_xi1xi1 +
      ( eig_xi1$vectors[,i] %*% t(eig_xi1$vectors[,i]) )
  }

  projector_xi1xi1
}

##############################################

.compute_GARCH_error_op <- function(X2, alpha, beta, delta, C_eps, C_X){
  error <- rep(0, dim(X2)[3])
  estimate <- .predict_GARCH_op(X2, alpha, beta, delta, C_eps)

  for(i in (2):dim(X2)[3]){
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
  for(i in (1+p):dims[3]){
    Sigma[,,i] <- delta +
      .apply_tensor(op = alpha,
                    X_apply = X2[,,i-j]) +
      .apply_tensor(op = beta,
                    X_apply = X2[,,i-j] %*% solve(C_eps))
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

