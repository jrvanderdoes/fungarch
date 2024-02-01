compare_methods_pointwise <- function(result_pw, result_or, X_est1){
  X_est <- X_est1
  resolution <- nrow(X_est)
  times <- seq(0,1,length.out=resolution)
  N <- ncol(X_est)

  # C_eps
  C_eps <- matrix(NA,resolution,resolution)
  times <- 1:resolution/resolution
  for (i in 1:resolution){
    C_eps[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
  }

  ## Get True X2 Values
  X2 <- X_est^2

  ## result_pwor(acle)
  estimate_pwor <- pointwise_next_step_ests(X2,result_or$hat_alpha,result_or$hat_delta)

  ## result_pw
  estimate_pw <- pointwise_next_step_ests(X2,result_pw$hat_alpha,result_pw$hat_delta)

  ## result_op
  # estimate_op <- .predictARCH_tensor_new(X2,
  #                                     alpha = result_op$hat_alpha,
  #                                     delta = result_op$hat_delta,
  #                                     C_eps = C_eps)

  ## Get Errors
  error <- data.frame('OR'=rep(0, ncol(X_est)),
                       'PW'=0,
                       'Mean'=0,
                       'Prev'=0)
  for(i in 2:ncol(X_est)){
    # TODO: Check integration
    error[i,] <- c(
      dot_integrate(abs(X2[,i] - estimate_pwor[,i])),
      dot_integrate(abs(X2[,i] - estimate_pw[,i])),
      dot_integrate(abs(X2[,i] - rowMeans(X2))),
      dot_integrate(abs(X2[,i] - X2[,i-1]))
    )
  }

  colSums(error)
}




get_oracle <- function(dat){
  s <- seq(0,1,length.out=nrow(dat$garch_mat))
  alpha_true <- matrix(nrow=nrow(dat$garch_mat),ncol=nrow(dat$garch_mat))
  for(s1 in 1:length(s)){
    for(s2 in 1:length(s)){
      alpha_true[s1,s2] <- dat$params$alpha_par(s[s1],s[s2])
    }
  }

  d <- seq(0,1,length.out=nrow(dat$garch_mat))
  delta_true <- dat$params$delta_par(d)

  list('hat_delta'=delta_true,
       'hat_alpha'=alpha_true)
}

##############################


pointwise_next_step_ests <- function(X2, alpha, delta){
  estimate <- matrix(NA, ncol = ncol(X2),nrow = nrow(X2))

  for(i in 2:ncol(X2)){
    # TODO:: Change to matrix mult
    for(j in 1:nrow(X2)){
      vector_op <- alpha[j, ] * X2[,(i-1)] # X^2
      estimate[j,i] <- delta[j] + dot_integrate(as.matrix(vector_op))
    }
  }

  estimate
}

operator_or_next_step_ests <- function(dat, alpha_par=NULL, delta_par=NULL){
  resolution <- nrow(dat)
  times <- seq(0,1,length.out=resolution)
  N <- ncol(dat)


  C_eps <- matrix(NA,resolution,resolution)
  for (i in 1:resolution){
    C_eps[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
  }

  if(is.null(alpha_par))
    alpha_par <- function(s,t,u,v){ s*(1-s)*t*(1-t)*u*(1-u)*v*(1-v) }
  if(is.null(delta_par))
    delta_par <- function(s,t){s*(1-s)*t*(1-t)}

  # Define variance and GARCH matrices
  estimate <- array(NA, c(resolution, resolution, N))

  vector_op1 <- vector_op <-
    array(dim = c(resolution,resolution,resolution))
  for(i in 2:N){
    prev_tensor <- dat[,i-1] %*% t(dat[,i-1])
    for(j in 1:resolution){
      # Over all t
      for(k in 1:resolution){
        vector_op[,,k] <-
          as.matrix(
            sapply(times, alpha_par,
                   s=times[j],t=times[k],u=times)
          ) %*% prev_tensor
      }
      # TODO:: Convert to trap int
      vector_op1[,,j] <- rowMeans(vector_op, dims=2)
    }
    # Formula (2)
    #  delta + sum alpha(X2) + sum beta(Sigma2)
    # TODO:: Convert to trap int
    sim_sigma2_matrix <-
      sapply(times,delta_par,t=times) + rowMeans(vector_op1, dims=2)

    #    See (66)
    Sigma2_svd <- svd(sim_sigma2_matrix)
    Sigma2_sqrt <- Sigma2_svd$u %*% diag(sqrt(Sigma2_svd$d)) %*% t(Sigma2_svd$v)
    # Formula (2)
    estimate[,,i] <- Sigma2_sqrt %*% C_eps %*% Sigma2_sqrt
  }

  estimate
}

operator_est_next_step_ests <- function(X2, alpha, delta){
  dims <- dim(X2)
  resolution <- dims[1]
  times <- seq(0,1,length.out=resolution)
  N <- dims[3]

  C_eps <- matrix(NA,resolution,resolution)
  for (i in 1:resolution){
    C_eps[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
  }


  Sigma <- array(NA, dims)

  for(i in 2:dims[3]){
    for(j in 1:dims[1]){
      vector_op <- alpha[j, ] * X2[,,(i-1)] # X^2
      Sigma[j,,i] <- delta[j,] + colSums(vector_op)/dims[1]
    }
  }

  Estimate <- array(NA, dims)
  for(i in 2:dims[3]){
    Sigma_svd <- svd(Sigma[,,i])
    Sigma_sqrt <- Sigma_svd$u %*% diag(sqrt(Sigma_svd$d)) %*% t(Sigma_svd$v)

    Estimate[,,i] <- Sigma_sqrt %*% (C_eps) %*% Sigma_sqrt
  }

  Estimate
}

compare_methods_operator <- function(result_pw, result_op, result_opor, X_est1){
  X_est <- X_est1
  resolution <- nrow(X_est)
  times <- seq(0,1,length.out=resolution)
  N <- ncol(X_est)

  # C_eps
  C_eps <- matrix(NA,resolution,resolution)
  times <- 1:resolution/resolution
  for (i in 1:resolution){
    C_eps[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
  }

  ## Get True X2 Values
  X2 <- array(apply(X_est, 2, tcrossprod),
                 dim=c(nrow(X_est),nrow(X_est),ncol(X_est)))
  # X2_op1 <- array(apply(X_est, 2, tcrossprod),
  #                 dim=c(nrow(X_est),nrow(X_est),ncol(X_est))) %*% C_eps

  ## result_pwor(acle)
  # estimate_pwor_pw <- pointwise_next_step_ests(
  #   X_est^2,result_pwor$hat_alpha,result_pwor$hat_delta)
  #
  # estimate_pwor_op <- array(NA, dim(X2_op))
  # estimate_pwor_op[,,1] <- replicate(n = nrow(X_est),
  #                                  expr = result_pwor$hat_delta,
  #                                  simplify = "dataframe")
  # estimate_pwor_op1 <- estimate_pwor_op
  # for(i in 2:ncol(X_est)){
  #   estimate_pwor_op[,,i] <- tcrossprod(estimate_pwor_pw[,i],estimate_pwor_pw[,i])
  #   estimate_pwor_op1[,,i] <- tcrossprod(estimate_pwor_pw[,i],estimate_pwor_pw[,i]) %*% C_eps
  # }

  ## result_opor(acle)
  estimate_or <-
    operator_or_next_step_ests(X_est,result_opor$hat_alpha,result_opor$hat_delta)


  ## result_pw
  estimate_pw_pw <- pointwise_next_step_ests(X_est^2,result_pw$hat_alpha,result_pw$hat_delta)
  estimate_pw_op <- array(NA, dim(X2))
  estimate_pw_op[,,1] <- replicate(n = nrow(X_est),
                                expr = result_pw$hat_delta,
                                simplify = "dataframe")
  estimate_pw_op1 <- estimate_pw_op
  for(i in 2:ncol(X_est)){
    estimate_pw_op[,,i] <- tcrossprod(estimate_pw_pw[,i],estimate_pw_pw[,i]) #%*% C_eps
    estimate_pw_op[,,i] <- (pmax(0, estimate_pw_op[,,i]))^(1/4)
    #tmp1 <- sqrt(pmax(0,estimate_pw_pw[,i]))
    #estimate_pw_op[,,i] <- tcrossprod(tmp1,tmp1)
  }

  ## result_op
  estimate_op <- operator_est_next_step_ests(X2,
                                         alpha = result_op$hat_alpha,
                                         delta = result_op$hat_delta)

  ## Get Errors
  #estimate_pwor_op, estimate_pwor_op1
  #estimate_pwor_pw
  #estimate_pw_op, estimate_pw_op1
  #estimate_op
  error <- data.frame('OR'=rep(0,ncol(X_est)),
                      'OP'=0,
                      'PW'=0,
                      'Mean'=0,
                      'Prev'=0)
  for(i in 2:ncol(X_est)){
    error[i-1,] <- c(
      dot_integrate(abs(X2[,,i] - estimate_or[,,i])), # OR
      dot_integrate(abs(X2[,,i] - estimate_op[,,i])), # OP
      dot_integrate(abs(X2[,,i] - estimate_pw_op[,,i])), # PW
      dot_integrate(abs(X2[,,i] - rowMeans(X2,dims = 2))), # Mean
      dot_integrate(abs(X2[,,i] - X2[,,i-1])) # Prev
    )
  }

  colSums(error)
}
