#' Estimate diagonal functional GARCH parameters
#'
#' Estimate diagonal autoregressive coefficients and the innovation covariance
#' using automatic basis-dimension selection.
#'
#' @param X Numeric matrix of functional observations, with observations in
#'   columns.
#' @param C_eps Innovation covariance matrix.
#' @param basis Basis-function matrix.
#' @param eigenvalues Eigenvalues associated with `basis`.
#' @param w Approximation threshold used for dimension selection and
#'   regularization.
#' @param p Number of autoregressive lags.
#'
#' @return A list containing:
#'   describe{
#'     \item{`alpha`}{Estimated diagonal autoregressive coefficients.}
#'     \item{`delta`}{Estimated innovation covariance contribution.}
#'     \item{`K_Ce`}{Number of retained innovation-covariance basis functions.}
#'     \item{`K`}{Number of retained fitting basis functions.}
#'   }
#'
#' @noRd
.fit_diag_model_new <- function(X,
                            C_eps, basis, eigenvalues,
                            w=0.05, p=1){#w1=0.1, p=1){#, dimen=nrow(C_eps)){
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X))) {
    stop("X must be a finite numeric matrix", call. = FALSE)
  }
  if (!is.matrix(C_eps) || !is.numeric(C_eps) ||
      nrow(C_eps) != ncol(C_eps) || nrow(C_eps) != nrow(X) ||
      any(!is.finite(C_eps))) {
    stop("C_eps must be a finite square matrix matching X", call. = FALSE)
  }
  if (!is.matrix(basis) || !is.numeric(basis) ||
      nrow(basis) != nrow(X) || ncol(basis) < 1L ||
      any(!is.finite(basis))) {
    stop("basis must have matching rows and at least one finite numeric column", call. = FALSE)
  }
  if (!is.numeric(eigenvalues) || length(eigenvalues) < 1L ||
      any(!is.finite(eigenvalues)) || any(eigenvalues <= 0)) {
    stop("eigenvalues must contain at least one positive finite value", call. = FALSE)
  }
  if (length(w) != 1L || !is.numeric(w) || !is.finite(w) || w < 0 || w > 1) {
    stop("w must be a single finite value between 0 and 1", call. = FALSE)
  }
  if (length(p) != 1L || !is.numeric(p) || !is.finite(p) ||
      p < 1 || p != floor(p)) {
    stop("p must be a single positive integer", call. = FALSE)
  }
  p <- as.integer(p)
  if (ncol(X) <= p + 1L) {
    stop("X must contain more than p + 1 observations", call. = FALSE)
  }

  #################
  ## Vars
  #################
  N <- ncol(X)
  n <- nrow(X)
  w1 <- w


  #################
  ## Dimension Selection
  #################
  # Fall back to the full basis if the error threshold is never reached.
  max_basis <- min(n, ncol(basis))
  dimen <- max_basis
  Xd <- matrix(0,nrow=n ,ncol=N)
  for(j in seq_len(max_basis)){
    tmp <- tmp1 <- rep(0,N)
    for(i in 1:N){
      Xd[,i] <- Xd[,i] + t((X[,i] %*% basis[,j]) %*% basis[,j]) / n

      tmp[i] <- sum((X[,i] - Xd[,i])^2)
      tmp1[i] <- sum((X[,i])^2)
    }

    if(sum(tmp)/sum(tmp1) < w1) {
      dimen = j
      break
    }
  }

  if (length(eigenvalues) < dimen) {
    stop("eigenvalues must include one value for each selected basis function", call. = FALSE)
  }

  #################
  ## Compute Matrices
  #################
  ## C_eps and related
  K_Ce <- dimen#which.min(cumsum(eigenvalues)/sum(eigenvalues) < 1-w)

  C_eps_dag <- matrix(0, nrow=nrow(C_eps), ncol=ncol(C_eps))
  for(i in 1:K_Ce){
    C_eps_dag <- C_eps_dag +
      ( 1 / eigenvalues[i] ) * (basis[,i] %o% basis[,i] )
  }
  # Note C_eps %*% C_eps_dag for max K_Ce is n^2*Identity




  tmp <- .CD_hat_new(X, p, dimen, C_eps_dag, basis)

  alpha_tmp <- tryCatch({
    tmp$D_hat %*% solve(tmp$C_hat)
  }, error = function(e){
    # warning('C_hat not invertible')
    eigs <- eigen(tmp$C_hat)
    dimen_inv <- which.min(cumsum(eigs$values)/sum(eigs$values) < 1 - w)
    C_hat_dag <- matrix(0, nrow=nrow(eigs$vectors), ncol = nrow(eigs$vectors))
    for(i in 1:dimen_inv){
      C_hat_dag <- C_hat_dag +
        ( 1 / eigs$values[i] ) * (eigs$vectors[,i] %o% eigs$vectors[,i] )
    }

    tmp$D_hat %*% C_hat_dag
  })

  if(dimen>1){
    alpha <-
      sapply(1:p, function(idx){
        diag(alpha_tmp[,(idx-1)*dimen + 1:dimen])
      })
  }else{
    alpha <- alpha_tmp
  }

  X2 <- array(apply(X, 2, tcrossprod),
              dim=c(nrow(X),nrow(X),ncol(X)))
  m_X2 <- rowMeans(X2, dims = 2)

  tmp_alpha <- sapply(1:dimen, function(d){
    int_vals <- (t( m_X2 %*% basis[,d] ) %*% basis[,d])[[1]] / n^2
    phi_front <- ( basis[,d] %*% t(basis[,d]) )

    # Compute then add all alphas parts (e.g. over p) at given dimension
    tmp_sums <- sapply(1:p, function(idx){
      alpha[d, idx] * phi_front * int_vals
    },simplify = F)

    Reduce('+', tmp_sums)
  },simplify = F)
  alpha_mp <- Reduce('+', tmp_alpha)

  list('alpha'=alpha,
       'delta'= matrix(
         pmax(0, (m_X2 - (alpha_mp %*% C_eps)/n) %*% C_eps_dag / n),
         nrow=nrow(X)
        ),
       'K_Ce'=K_Ce,
       'K' = dimen)
}

#################
## C/D
#################
#' Estimate covariance and cross-covariance matrices for the new diagonal fit.
#'
#' @param X Numeric matrix of observations.
#' @param p Number of autoregressive lags.
#' @param dimen Number of retained basis functions.
#' @param C_eps_dag Truncated inverse innovation covariance.
#' @param basis Basis-function matrix.
#'
#' @return A list containing `C_hat`, the lagged predictor covariance matrix,
#'   and `D_hat`, the corresponding cross-covariance matrix.
#'
#' @noRd
.CD_hat_new <- function(X, p, dimen, C_eps_dag, basis){
  N <- ncol(X)
  r <- nrow(X)
  Chat <- Dhat <- matrix(0, nrow=dimen, ncol=dimen)

  #############
  # Common
  # B <- C_eps_dag #%*% basis  # Tried all transposed of basis and div by r

  #############
  # METHOD
  M <- array(dim=c(dimen,dimen,N))
  for(i in 1:N){
    for(j in 1:dimen){
      for(k in 1:dimen){
        M[j,k,i] <- (X[,i] %*% basis[,j]) * (X[,i] %*% basis[,k]) / r^2
      }
    }
  }

  C_eps_dagM = matrix(NA,dimen,dimen)
  for(ii in 1:dimen){
    for(jj in 1:dimen){
      C_eps_dagM[ii,jj] <-
        fChange:::dot_integrate( C_eps_dag * ( basis[,ii] %*% t(basis[,jj]) ) )
    }
  }

  M_mean <- rowMeans(M,dim=2)
  for(i in 1:N){
    M[,,i] <- M[,,i] - M_mean
  }

  X_M <- matrix(,nrow=dimen, ncol=N)
  X_Mp <- matrix(,nrow=dimen*p, ncol=N)
  for(i in 1:N){
    ## TODO:: Fix for 1D
    # X^X2,[p]
    # print(dim(M))
    if(dimen>1){
      X_Mp[,i] <- #M[i]
        as.vector( sapply(1:p, function(idx){
          if(i-idx+1<=0) return(rep(NA,nrow(M)))
          diag(M[,,i-idx+1])
        }) )

      # X^X2_eps
      X_M[,i] <- diag(M[,,i] %*% C_eps_dagM)
    }else{
      X_Mp[,i] <- M[i]

      # X^X2_eps
      X_M[,i] <- diag(M[i] %*% C_eps_dagM)
    }

  }

  C_hat <- X_Mp[,N] %*% t(X_Mp[,N]) # Original p
  D_hat <- matrix(0,nrow=dimen,ncol=p*dimen)
  for(i in p:(N-1)){
    C_hat <- C_hat + X_Mp[,i] %*% t(X_Mp[,i])
    D_hat <- D_hat + X_M[,i+1] %*% t(X_Mp[,i])
  }
  C_hat <- C_hat / (N - p)
  D_hat <- D_hat / (N - p - 1)

  list('C_hat'= C_hat,
       'D_hat'= D_hat
  )
}

#' Predict the conditional covariance using the new diagonal fit.
#'
#' @param resolution Functional resolution.
#' @param alphas Estimated diagonal autoregressive coefficients.
#' @param basis Basis-function matrix.
#' @param delta Estimated innovation covariance.
#' @param data Recent functional observations.
#' @param C_eps Innovation covariance matrix.
#'
#' @return A predicted conditional covariance matrix.
#'
#' @noRd
.pred_diag_new <- function(resolution, alphas, basis, delta, data, C_eps){
  # OLD: The following duplicated alpha and covariance construction is
  # retained for comparison; it is now handled by .build_diag_sigma().
  # OLD: Each lag received the sum of all lag coefficients. This is retained
  # # for comparison while the lag-specific implementation is evaluated.
  # alpha_applied <- replicate(ncol(alphas),
  #                            matrix(0, nrow=resolution^2, ncol = resolution^2),
  #                            simplify=FALSE)
  # for(p in 1:ncol(alphas)){
  #   for(k in 1:nrow(alphas)){
  #     hkk <- as.vector(basis[,k] %o% basis[,k])
  #     alphas_list <- sapply(1:ncol(alphas), function(d){
  #       alphas[k,d] * (hkk %o% hkk)
  #     },simplify = F)
  #     alpha_applied[[p]] <- alpha_applied[[p]] + Reduce('+', alphas_list)
  #   }
  # }

  # alpha_applied <- lapply(seq_len(ncol(alphas)), function(p) {
  #   result <- matrix(0, nrow = resolution^2, ncol = resolution^2)
  #   for (k in seq_len(nrow(alphas))) {
  #     hkk <- as.vector(basis[, k] %o% basis[, k])
  #     result <- result + alphas[k, p] * (hkk %o% hkk)
  #   }
  #   result
  # })
  # sigma <- as.vector(delta) +
  #   Reduce('+', sapply(1:ncol(alphas), function(idx){
  #     operator(alpha_applied[[idx]],
  #              tcrossprod(data[, ncol(data)-idx+1],
  #                         data[, ncol(data)-idx+1]))
  #   }, simplify = FALSE))

  sigma <- .build_diag_sigma(resolution, alphas, basis, delta, data)

  # Get Sqrt(Sigma)
  sigma_svd <- svd(sigma)
  sigma_sqrt <- sqrt(resolution) *
    sigma_svd$u %*% diag(sqrt(sigma_svd$d)) %*% t(sigma_svd$v)

  # FUTURE ALTERNATIVE: Symmetrize the covariance before taking its square root.
  # sigma <- (sigma + t(sigma)) / 2
  # sigma_eig <- eigen(sigma, symmetric = TRUE)
  # tol <- .Machine$double.eps * max(1, max(abs(sigma_eig$values)))
  # if (any(sigma_eig$values < -tol)) {
  #   stop("conditional covariance is not positive semidefinite", call. = FALSE)
  # }
  # sigma_eig$values <- pmax(sigma_eig$values, 0)
  # sigma_sqrt <- sqrt(resolution) *
  #   sigma_eig$vectors %*% diag(sqrt(sigma_eig$values)) %*%
  #   t(sigma_eig$vectors)

  # Predict
  sigma_sqrt %*% C_eps %*% sigma_sqrt / resolution^2
}

#' Compute the square root of the conditional covariance for the new fit.
#'
#' @inheritParams .pred_diag_new
#'
#' @return A matrix square root of the conditional covariance.
#'
#' @noRd
.sigma2_diag <- function(resolution, alphas, basis, delta, data, C_eps){
  # OLD: The following duplicated alpha and covariance construction is
  # retained for comparison; it is now handled by .build_diag_sigma().
  # OLD: Each lag received the sum of all lag coefficients. This is retained
  # # for comparison while the lag-specific implementation is evaluated.
  # alpha_applied <- replicate(ncol(alphas),
  #                            matrix(0, nrow=resolution^2, ncol = resolution^2),
  #                            simplify=FALSE)
  # for(p in 1:ncol(alphas)){
  #   for(k in 1:nrow(alphas)){
  #     hkk <- as.vector(basis[,k] %o% basis[,k])
  #     alphas_list <- sapply(1:ncol(alphas), function(d){
  #       alphas[k,d] * (hkk %o% hkk)
  #     },simplify = F)
  #     alpha_applied[[p]] <- alpha_applied[[p]] + Reduce('+', alphas_list)
  #   }
  # }

  # alpha_applied <- lapply(seq_len(ncol(alphas)), function(p) {
  #   result <- matrix(0, nrow = resolution^2, ncol = resolution^2)
  #   for (k in seq_len(nrow(alphas))) {
  #     hkk <- as.vector(basis[, k] %o% basis[, k])
  #     result <- result + alphas[k, p] * (hkk %o% hkk)
  #   }
  #   result
  # })
  # sigma <- as.vector(delta) +
  #   Reduce('+', sapply(1:ncol(alphas), function(idx){
  #     operator(alpha_applied[[idx]],
  #              tcrossprod(data[, ncol(data)-idx+1],
  #                         data[, ncol(data)-idx+1]))
  #   }, simplify = FALSE))

  sigma <- .build_diag_sigma(resolution, alphas, basis, delta, data)

  # Get Sqrt(Sigma)
  sigma_svd <- svd(sigma)
  sigma_sqrt <- sqrt(resolution) *
    sigma_svd$u %*% diag(sqrt(sigma_svd$d)) %*% t(sigma_svd$v)

  # FUTURE ALTERNATIVE: Symmetrize the covariance before taking its square root.
  # sigma <- (sigma + t(sigma)) / 2
  # sigma_eig <- eigen(sigma, symmetric = TRUE)
  # tol <- .Machine$double.eps * max(1, max(abs(sigma_eig$values)))
  # if (any(sigma_eig$values < -tol)) {
  #   stop("conditional covariance is not positive semidefinite", call. = FALSE)
  # }
  # sigma_eig$values <- pmax(sigma_eig$values, 0)
  # sigma_sqrt <- sqrt(resolution) *
  #   sigma_eig$vectors %*% diag(sqrt(sigma_eig$values)) %*%
  #   t(sigma_eig$vectors)

  sigma_sqrt
}

#' Build the conditional covariance used by the diagonal prediction helpers.
#'
#' @param resolution Functional resolution.
#' @param alphas Estimated diagonal autoregressive coefficients.
#' @param basis Basis-function matrix.
#' @param delta Estimated innovation covariance.
#' @param data Recent functional observations.
#'
#' @return A conditional covariance matrix before applying `C_eps`.
#'
#' @noRd
.build_diag_sigma <- function(resolution, alphas, basis, delta, data){
  alpha_applied <- lapply(seq_len(ncol(alphas)), function(p) {
    result <- matrix(0, nrow = resolution^2, ncol = resolution^2)
    for (k in seq_len(nrow(alphas))) {
      hkk <- as.vector(basis[, k] %o% basis[, k])
      result <- result + alphas[k, p] * (hkk %o% hkk)
    }
    result
  })

  sigma <- as.vector(delta) +
    Reduce('+', lapply(seq_len(ncol(alphas)), function(idx) {
      operator(
        alpha_applied[[idx]],
        tcrossprod(data[, ncol(data) - idx + 1],
                   data[, ncol(data) - idx + 1])
      )
    }))

  matrix(sigma, nrow = resolution, ncol = resolution)
}



#############
#' Estimate diagonal functional GARCH parameters with Tikhonov regularization.
#'
#' This is the regularized variant of `.fit_diag_model_new()` used by the
#' experimental workflows.
#'
#' @param X Numeric matrix of functional observations.
#' @param C_eps Innovation covariance matrix.
#' @param basis Basis-function matrix.
#' @param eigenvalues Eigenvalues associated with `basis`.
#' @param w Regularization and dimension-selection threshold.
#' @param p Number of autoregressive lags.
#' @param lambda Optional starting Tikhonov regularization parameter. If `NA`,
#'   a default based on the leading eigenvalue is used.
#' @param training Fraction of observations used for fitting. Values below 1
#'   reserve the remaining observations for prediction-error tuning.
#' @param alp Quantile level used when tuning the regularization parameter.
#'
#' @return A list containing `alpha`, the fitted diagonal coefficients,
#'   `delta`, the estimated innovation covariance contribution, and `K_Ce`,
#'   the number of retained innovation-covariance basis functions.
#'
#' @noRd
.fit_diag_model_tik <- function(X,
                                C_eps, basis, eigenvalues,
                                w=0.05, p=1, #dimen=nrow(C_eps), #w1=0.1,
                                lambda=NA, training=0.8, alp=0.05){
  if (length(training) != 1L || !is.numeric(training) ||
      !is.finite(training) || training <= 0) {
    stop("training must be a single positive finite value", call. = FALSE)
  }
  w1 <- w

  #################
  ## Dimension Selection
  #################
  max_basis <- min(nrow(X), ncol(basis))
  dimen <- max_basis
  Xd <- matrix(0,nrow=nrow(X), ncol=ncol(X))
  for(j in seq_len(max_basis)){
    tmp <- rep(0,ncol(X))
    for(i in 1:ncol(X)){
      Xd[,i] <- Xd[,i] + t((X[,i] %*% basis[,j]) %*% basis[,j]) / nrow(X)

      tmp[i] <- sum(abs(X[,i] - Xd[,i])) / sum(abs(X[,i]))
    }

    if(mean(tmp) < w1) {
      dimen = j
      break
    }
  }

  #################
  ## Vars
  #################
  # Inv C_eps
  ## C_eps and related
  K_Ce <- dimen#which.min(cumsum(eigenvalues)/sum(eigenvalues) < 1-w)

  C_eps_dag <- matrix(0, nrow=nrow(C_eps), ncol=ncol(C_eps))
  for(i in 1:K_Ce){
    C_eps_dag <- C_eps_dag +
      ( 1 / eigenvalues[i] ) * (basis[,i] %o% basis[,i] )
  }

  ## USE PART OF DATA FOR FIT
  if(training < 1){
    # OLD:
    # X_part <- X[,1:(training * ncol(X))]
    # X_pred <- X[,(training * ncol(X)+1):ncol(X)]
    n_train <- floor(training * ncol(X))
    if (n_train <= p + 1L || n_train >= ncol(X)) {
      stop("training must leave enough observations for fitting and prediction", call. = FALSE)
    }
    X_part <- X[, seq_len(n_train), drop = FALSE]
    X_pred <- X[, seq.int(n_train + 1L, ncol(X)), drop = FALSE]
  }else{
    X_part <- X
    X_pred <- X
  }
  N <- ncol(X_part)
  n <- nrow(X)

  # #################
  # ## C/D
  # #################
  # .CD_hat <- function(X, p, dimen, C_eps_dag, basis){
  #   N <- ncol(X)
  #   r <- nrow(X)
  #   Chat <- Dhat <- matrix(0, nrow=dimen, ncol=dimen)
  #
  #   #############
  #   # METHOD
  #   M <- array(dim=c(dimen,dimen,N))
  #   for(i in 1:N){
  #     for(j in 1:dimen){
  #       for(k in 1:dimen){
  #         M[j,k,i] <- (X[,i] %*% basis[,j]) * (X[,i] %*% basis[,k]) / r^2
  #       }
  #     }
  #   }
  #
  #   C_eps_dagM = matrix(NA,dimen,dimen)
  #   for(ii in 1:dimen){
  #     for(jj in 1:dimen){
  #       C_eps_dagM[ii,jj] <-
  #         fChange:::dot_integrate( C_eps_dag * ( basis[,ii] %*% t(basis[,jj]) ) )
  #     }
  #   }
  #
  #   M_mean <- rowMeans(M,dim=2)
  #   for(i in 1:N){
  #     M[,,i] <- M[,,i] - M_mean
  #   }
  #
  #   X_M <- matrix(,nrow=dimen, ncol=N)
  #   X_Mp <- matrix(,nrow=dimen*p, ncol=N)
  #   for(i in 1:N){
  #     ## TODO:: Fix for 1D
  #     # X^X2,[p]
  #     X_Mp[,i] <- #M[i]
  #       as.vector( sapply(1:p, function(idx){
  #         if(i-idx+1<=0) return(rep(NA,nrow(M)))
  #         diag(M[,,i-idx+1])
  #       }) )
  #
  #     # X^X2_eps
  #     X_M[,i] <- diag(M[,,i] %*% C_eps_dagM)
  #   }
  #
  #   C_hat <- X_Mp[,N] %*% t(X_Mp[,N]) # Original p
  #   D_hat <- matrix(0,nrow=dimen,ncol=p*dimen)
  #   for(i in p:(N-1)){
  #     C_hat <- C_hat + X_Mp[,i] %*% t(X_Mp[,i])
  #     D_hat <- D_hat + X_M[,i+1] %*% t(X_Mp[,i])
  #   }
  #   C_hat <- C_hat / (N - p)
  #   D_hat <- D_hat / (N - p - 1)
  #
  #   list('C_hat'= C_hat,
  #        'D_hat'= D_hat
  #   )
  # }

  #################
  ## OPTIM
  #################
  lambda[is.na(lambda)] <- 1/2 * eigenvalues[1]

  X2_part <- array(apply(X_part, 2, tcrossprod),
                   dim=c(nrow(X_part),nrow(X_part),ncol(X_part)))
  m_X2_part <- rowMeans(X2_part, dims = 2)

  optim_results <-
    optim(lambda, function(lambda, X1, N, n, projector,
                              m_X2_1, basis, X2_1, alp,eigenvalues){

      ## C, D
      tmp <- .CD_hat_new(X1, p, dimen, C_eps_dag, basis)

      ## Alpha
      alpha_tmp <- tmp$D_hat %*%
        solve(tmp$C_hat + diag(lambda, nrow(tmp$C_hat), nrow(tmp$C_hat)) ) #%*% projector

      if(dimen>1){
        alpha <-
          sapply(1:p, function(idx){
            diag(alpha_tmp[,(idx-1)*dimen + 1:dimen])
          })
      }else{
        alpha <- alpha_tmp
      }

      tmp_alpha <- sapply(1:dimen, function(d){
        int_vals <- (t( m_X2_1 %*% basis[,d] ) %*% basis[,d])[[1]] / n^2
        phi_front <- ( basis[,d] %*% t(basis[,d]) )

        # Compute then add all alphas parts (e.g. over p) at given dimension
        tmp_sums <- sapply(1:p, function(idx){
          alpha[d, idx] * phi_front * int_vals ## TODO:: REVIEW d was dimen
        },simplify = F)

        Reduce('+', tmp_sums)
      },simplify = F)
      alpha_mp <- Reduce('+', tmp_alpha)

      ## Delta
      delta <- matrix(
        pmax(0, (m_X2_part - (alpha_mp %*% C_eps)/n) %*% C_eps_dag / n),
        nrow=nrow(X1)
      )

      ## Errors (check Func)
      sds <- apply(X1,MARGIN = 1,sd)

      err <- rep(0,ncol(X_pred))
      estimate <- .pred_diag_new(n, alpha, basis, delta,
                                 X1, C_eps)
      quant <- sqrt(pmax(0,diag(estimate))) * qnorm(alp,sd=sds)
      u <- X2_1[,1] - quant
      err[1] <- fChange:::dot_integrate(ifelse(u>=0,alp*abs(u),(1-alp)*abs(u)))

      for(i in 2:(ncol(X_pred)-1)){
        estimate <- .pred_diag_new(n, alpha, basis, delta,
                                   cbind(X1, X2_1[,1:(i-1)]), C_eps)
        quant <- sqrt(pmax(0,diag(estimate))) * qnorm(alp, sd=sds)
        u <- X2_1[,i] - quant
        err[i] <- fChange:::dot_integrate(ifelse(u>=0,alp*abs(u),(1-alp)*abs(u)))
      }

      sum(err)
    },
    method="L-BFGS-B",
    lower=10^-7,upper=500000, # change upper
    X1=X_part, N=N, n=n,m_X2_1=m_X2_part,
    basis=basis, eigenvalues=eigenvalues,
    X2_1=X_pred, alp=alp)

    lambda <- optim_results$par

    ### FINAL
    tmp <- .CD_hat_new(X, p, dimen, C_eps_dag, basis)
    alpha_tmp <- tmp$D_hat %*%
      solve(tmp$C_hat + diag(lambda, nrow(tmp$C_hat), nrow(tmp$C_hat)) ) #%*% projector

    ## alpha
    if(dimen>1){
      alpha <-
        sapply(1:p, function(idx){
          diag(alpha_tmp[,(idx-1)*dimen + 1:dimen])
        })
    }else{
      alpha <- alpha_tmp
    }


    X2 <- array(apply(X, 2, tcrossprod),
                dim=c(nrow(X),nrow(X),ncol(X)))
    m_X2 <- rowMeans(X2, dims = 2)

    tmp_alpha <- sapply(1:dimen, function(d){
      int_vals <- (t( m_X2 %*% basis[,d] ) %*% basis[,d])[[1]] / n^2
      phi_front <- ( basis[,d] %*% t(basis[,d]) )

      # Compute then add all alphas parts (e.g. over p) at given dimension
      tmp_sums <- sapply(1:p, function(idx){
        alpha[d, idx] * phi_front * int_vals ## TODO:: REVIEW d was dimen
      },simplify = F)

      Reduce('+', tmp_sums)
    },simplify = F)
    alpha_mp <- Reduce('+', tmp_alpha)

    ## delta
    delta <- matrix(
      pmax(0, (m_X2 - (alpha_mp %*% C_eps)/n) %*% C_eps_dag / n),
      nrow=nrow(X)
    )

    list('alpha'=alpha,
         'delta'= delta,
         'K_Ce'=K_Ce)
}
