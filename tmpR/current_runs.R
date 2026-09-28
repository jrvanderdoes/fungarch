data <- generateFunctionalCondVol(21,500,'arch')
res <- operatorARCH(X=data$garch_mat,p=1)
or <- .get_oracle(data)

# TODO:: Discuss differences
data <- generateFunctionalCondVol_ff(21,1000,'arch')
res <- operatorARCH(X=data$garch_mat,p=1)
or <- .get_oracle_op(data)

plot_fd(res$hat_alpha)
plot_fd(or$hat_alpha)
plot_fd(res$hat_delta)
plot_fd(or$hat_delta)

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
operatorARCH <- function(X, w=0.001, thetas=NA, training=0.8, p=1){
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
    X2_part_eps[,,i] <- X2_part[,,i] %*% C_eps
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

  C_hat <- tmp$C_hat # was t()
  D_hat <- tmp$D_hat

  # After Eq 36 and Eq 39
  eig_CC <- eigen( C_hat %*% t(C_hat) )
  eig_Ce <- eigen( C_eps )
  K_CC <- which.min(cumsum(eig_CC$values)/sum(eig_CC$values) < 1-w)
  K_Ce <- which.min(cumsum(eig_Ce$values)/sum(eig_Ce$values) < 1-w)
  # Adjust since these are finite and not functions. Need integral to be 1
  eig_CC$vectors <- eig_CC$vectors * sqrt(n)
  eig_Ce$vectors <- eig_Ce$vectors * sqrt(n)

  # Compute Projections
  projector_CC <- matrix(0, nrow=n*p, ncol=n*p)
  projector_Ce <- matrix(0, nrow=n, ncol=n)
  for(i in 1:K_CC){
    projector_CC <- projector_CC +
      ( eig_CC$vectors[,i] %*% t(eig_CC$vectors[,i]) )
  }
  for(i in 1:K_Ce){
    projector_Ce <- projector_Ce +
      ( eig_Ce$vectors[,i] %*% t(eig_Ce$vectors[,i]) )
  }

  # Find the optimal Tikhonov parameters
  if(length(thetas)!=2 ||
     sum(is.na(thetas)) > 0 ){

    optim_results <-
      optim(c(10,10), fn=function(theta0, C_hat, D_hat,
                                  projector_CC, projector_Ce,
                                  n, p1, X2, C_X1){
        theta <- theta0[1]
        vu <- theta0[2]

        ests <- .estimate_params_ARCH_op(theta = theta, vu=vu, C_hat = C_hat,
                                         D_hat = D_hat, C_eps=C_eps, C_X=C_X1,
                                         projector_CC = projector_CC,
                                         projector_Ce = projector_Ce,
                                         n = n, p = p1)

        .compute_ARCH_error_op(
          X2 = X2, alpha = ests$alpha_hat, delta = ests$delta_hat,
          C_eps=C_eps, p=p1, C_X=C_X1)
      },
      method="L-BFGS-B",#"Brent",
      lower=10^-8,upper=500000000,
      C_hat=C_hat, D_hat=D_hat, projector_CC=projector_CC,
      projector_Ce=projector_Ce, n=n, p1=p, X2=X2_pred, C_X1=C_part_X)

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
    X2_eps[,,i] <- X2[,,i] %*% C_eps
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

  C_hat <- (tmp$C_hat) # was t()
  D_hat <- (tmp$D_hat)

  eig_CC <- eigen( C_hat %*% t(C_hat) )
  eig_Ce <- eigen( C_eps )
  K_CC <- which.min(cumsum(eig_CC$values)/sum(eig_CC$values) < 1-w)
  K_Ce <- which.min(cumsum(eig_Ce$values)/sum(eig_Ce$values) < 1-w)
  eig_CC$vectors <- eig_CC$vectors * sqrt(n)
  eig_Ce$vectors <- eig_Ce$vectors * sqrt(n)

  # Compute Projection
  projector_CC <- matrix(0, nrow=n*p, ncol=n*p)
  projector_Ce <- matrix(0, nrow=n, ncol=n)
  for(i in 1:K_CC){
    projector_CC <- projector_CC +
      ( eig_CC$vectors[,i] %*% t(eig_CC$vectors[,i]) )
  }
  for(i in 1:K_Ce){
    projector_Ce <- projector_Ce +
      ( eig_Ce$vectors[,i] %*% t(eig_Ce$vectors[,i]) )
  }

  ests <- .estimate_params_ARCH_op(theta = thetas[1], vu=thetas[2],
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


.estimate_params_ARCH_op <- function(theta, vu, C_hat,D_hat, C_eps, C_X,
                                     projector_CC,projector_Ce,
                                     n,p){
  # TODO:: Confirm t()

  # Compute alpha, Eq 36
  C_hat_st <- t(Conj(C_hat))
  C_dag <- C_hat_st %*% solve( C_hat %*% C_hat_st + diag(theta,n*p) )
  alpha_hat <- D_hat %*% C_dag %*% projector_CC
  alpha_hat  <- apply(alpha_hat, 1, pmax, 0)

  ## Delta, Eq 39
  C_eps_diag <- solve( C_eps + diag(vu, n) )
  delta_hat <-
    (C_X - .apply_tensor(alpha_hat,
                         t(do.call(cbind, replicate(p, C_X, simplify=FALSE)))
    ) %*%
      C_eps ) %*%
    C_eps_diag  %*%
    projector_Ce
  delta_hat <- apply(delta_hat, 1, pmax, 0)

  list('alpha_hat'=alpha_hat,
       'delta_hat'=delta_hat)
}

.estim_CD_hat <- function(X, Xp, Xpe, p){
  if(prod(dim(Xp) != dim(Xpe)))
    stop('Error: Xp and Xpe have different dimensions')

  N <- dim(Xp)[3]

  Chat <- matrix(0, nrow=nrow(Xp),ncol=nrow(Xp))
  Dhat <- matrix(0, nrow=ncol(Xp),ncol=nrow(Xp))

  # After Eq 35
  Cxp <- .mean_C(Xp,p=p,lag=0,lag2=0)
  Cxpe <- .mean_C(Xpe,p=p,lag=0,lag2=0)
  Cxppr <- .mean_C(Xp,p=p,lag=1,lag2=0)
  Cxpr <- .mean_C(X,p=p,lag=1,lag2=1)

  # C_hat, Eq 34
  for(k in p:N){
    Chat <- Chat +
      .tensor(Xp[,,k]-Cxp, Xpe[,,k]-Cxpe)
  }

  # D_hat, Eq 35
  for(k in p:(N-1)){
    Dhat <- Dhat +
      .tensor(Xp[,,k]-Cxppr, X[,,k+1]-Cxpr)
  }

  list('C_hat'= Chat / (N-p),
       'D_hat'= Dhat / (N-p-1)
  )
}


.tensor <- function(X,Y){
  Y %*% t(X) #/ nrow(Y)
}
.apply_tensor <- function(op, X){
  t(X) %*% op
}
# .apply_tensor <- function(op, X){
#   op %*% t(X) / ncol(op)
#   # X %*% t(op) / ncol(op)
#
# }

.mean_C <- function(X, p, lag,lag2=0){
  dims <- dim(X)
  mat <- matrix(0,nrow=dims[1],ncol=dims[2])

  for(k in p:(dims[3]-lag)){
    mat <- mat + X[,,k+lag2]
  }

  1/(dims[3]-p+1-lag) * mat
}

##############################################

.compute_ARCH_error_op <- function(X2, alpha, delta, C_eps, C_X, p){
  error <- rep(0, dim(X2)[3])
  estimate <- .predict_ARCH_op(X2,alpha,delta, C_eps, C_X, p)

  for(i in (1+p):dim(X2)[3]){
    # error[i] <- .approx_int(abs(X2[,,i] - estimate[,,i]))
    # error[i] <- .approx_int(rowSums(abs(X2[,,i] - estimate[,,i]))/ncol(X2))
    error[i] <- .approx_int(apply(abs(X2[,,i] - estimate[,,i]),
                                  MARGIN = 2,
                                  FUN = .approx_int))
  }

  sum(error)
}

.predict_ARCH_op <- function(X2, alpha, delta, C_eps, C_X, p){
  dims <- dim(X2)
  Sigma <- array(NA, dims)

  # Eq 2
  for(i in (1+p):dims[3]){
    Sigma[,,i] <- delta
    for(j in 1:p){
      Sigma[,,i] <- Sigma[,,i] +
        .apply_tensor(alpha[1:dims[1]+(j-1)*dims[1],], X2[,,i-j])
    }
  }

  # Eq 29
  Estimate <- array(NA, dims)
  for(i in (1+p):dims[3]){
    Sigma_svd <- svd(Sigma[,,i])
    Sigma_sqrt <- Sigma_svd$u %*% diag(sqrt(Sigma_svd$d)) %*% t(Sigma_svd$v)

    tmp <- matrix(nrow = dims[1], ncol = dims[2])
    for(j in 1:p){
      tmp <- tmp +
        .apply_tensor(alpha[1:dims[1]+(j-1)*dims[1],], X2[,,i-j]-C_X)
    }
    tmp <- tmp %*% C_eps # TODO:: Confirm??

    Estimate[,,i] <- Sigma_sqrt %*% (C_eps) %*% Sigma_sqrt + tmp
  }

  Sigma
}

