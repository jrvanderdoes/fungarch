generate_opARCH <- function(resolution=20, N=750,
                            alpha_par=function(s,t,u,v){ 1/2*s*(1-s)*t*(1-t)*u*(1-u)*v*(1-v) },
                            delta_par=function(s,t){s*(1-s)*t*(1-t)},
                            burnin=1000, p=5){
  # Define variables
  times <- seq(0,1,length.out=resolution)
  NBurn <- burnin + N
  # sim_sigma2_matrix <- array(NA, c(resolution, resolution, NBurn))
  sim_sigma2_matrix1 <- matrix(0,ncol = NBurn,nrow=resolution^2)
  arch_data <- arch_data1 <-  matrix(NA,nrow=resolution,ncol=NBurn)

  eps <- .generateFunctionalIID(resolution = resolution, N = NBurn)

  ## Define variance and GARCH matrices
  # Initialize
  for(i in 1:p){
    arch_data1[,i] <-
      sapply(times,delta_par, t=times) %*% eps[,i] / resolution
  }

  # # Attempt 1: Precompute Alpha parameters to use as needed
  # vector_op1 <- vector_op <-
  #   array(dim = c(resolution,resolution,resolution))
  # Try For speed
  vector_list <- list()
  for(j in 1:resolution){
    vector_list[[j]] <- list()
    # Over all t
    for(k in 1:resolution){
      # Makes first index of list s
      #  Makes second index of list t
      #  Makes rows of each u
      #  Makes cols of each v
      vector_list[[j]][[k]] <-
       as.matrix(
         sapply(times, alpha_par,
                s=times[j],t=times[k],u=times)
       )
    }
  }

  # Attempt 2: Compute Alpha matrix
  idx <- 0
  alpha <- matrix(nrow=resolution^2,ncol=resolution^2)
  for(j in 1:resolution){
    for(k in 1:resolution){
      idx <- idx + 1
      alpha[,idx] <- as.vector(t(vector_list[[j]][[k]]))
    }
  }

  # Set matrices
  # Sigma2_list <- Sigma2_sqrt_list <- list()
  Sigma2_list1 <- Sigma2_sqrt_list1 <- list()

  # Compute values
  for(i in (p+1):NBurn){
    # tmp_alpha <- array(0,dim = c(resolution,resolution))
    tmp_alpha1 <- rep(0, resolution^2)
   for(pp in 1:p){
     prev_tensor <-
       as.vector( arch_data1[,i-pp] %*% t(arch_data1[,i-pp]) )

     # Attempt 2:
     tmp_alpha1 <- tmp_alpha1 + alpha %*% prev_tensor

     # # Attempt 1:
     # for(j in 1:resolution){
     #   # Over all t
     #   for(k in 1:resolution){
     #     vector_op[,,k] <- vector_list[[j]][[k]] %*% prev_tensor / resolution
     #   }
     #   vector_op1[,,j] <- rowMeans(vector_op, dims=2)
     # }
     # # sum alpha(X2)
     # tmp_alpha <- tmp_alpha + rowMeans(vector_op1, dims=2)
   }

   # #  Attempt 1: delta + sum alpha(X2) + sum beta(Sigma2)
   # sim_sigma2_matrix[,,i] <-
   #   sapply(times,delta_par,t=times) +
   #   tmp_alpha
   #  Attempt 2: delta + sum alpha(X2) + sum beta(Sigma2)
   sim_sigma2_matrix1[,i] <-
     as.vector(sapply(times,delta_par,t=times)) +
     tmp_alpha1

   # # Attempt 1: Compute the sqrt matrix
   # Sigma2_svd <- svd(sim_sigma2_matrix[,,i])
   # Sigma2_sqrt <- sqrt(resolution) * Sigma2_svd$u %*%
   #   diag(sqrt(Sigma2_svd$d)) %*%
   #   t(Sigma2_svd$v)
   # arch_data[,i] <- Sigma2_sqrt %*% eps[,i] / resolution
   #
   # Sigma2_sqrt_list[[i]] <- Sigma2_sqrt
   # Sigma2_list[[i]] <- sim_sigma2_matrix[,,i]

   # Attempt 2: Compute the sqrt matrix
   tryCatch({
     Sigma2_svd1 <- svd(matrix(sim_sigma2_matrix1[,i],
                               ncol=resolution,nrow=resolution))
   }, error=function(e){
     View(sim_sigma2_matrix1[,i])
     View(sim_sigma2_matrix1)
     stop('Error')
     })
   Sigma2_sqrt1 <- sqrt(resolution) * Sigma2_svd1$u %*%
     diag(sqrt(Sigma2_svd1$d)) %*% t(Sigma2_svd1$v)
   arch_data1[,i] <- Sigma2_sqrt1 %*% eps[,i] / resolution
  }

  # View(arch_data)
  # View(arch_data1)

  list(data = arch_data1[, burnin + 1:N],
       sigma_mat = sim_sigma2_matrix1[, burnin + 1:N],
       error_matrix = eps[, burnin + 1:N],
       params=list(p=p,q=q,
                   alpha_par=alpha_par,
                   alpha = alpha,
                   delta_par=delta_par,
                   delta = as.vector(sapply(times,delta_par,t=times)),
                   Sigma2_list= Sigma2_list1[burnin + 1:N],
                   Sigma2_sqrt_list= Sigma2_sqrt_list1[burnin + 1:N]))#,
                   # vector_list=vector_list))
}


op5_alphaComparision <-
  function(
    resolution = 20,
    alpha = paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\',
                   'ARCH5send',3,'.rds'),
    delta = function(s,t){ pmin(t,s) },
    nSims=50,
    seed=NA){

    if(!is.na(seed)) set.seed(seed)

    p <- 2
    lengths <- c(50,100,250,500)
    error_results_yw <- error_results_yw1 <- error_results_yw2 <-
      error_results_mp <- error_results_mp1 <- error_results_mp2 <-
      data.frame('Mean'=rep(NA,length(lengths)), 'Var'=NA)
    errors_list_yw <- errors_list_mp <- list()

    for(i in 1:length(lengths)){
      cat(paste0('Length ',lengths[i],': '))
      errors_yw <- errors_yw1 <- errors_yw2 <-
        errors_mp <- errors_mp1 <- errors_mp2 <- rep(NA, nSims)
      for(sim in 1:nSims){
        cat(paste0(sim,', '))
        # dat <- generate_opARCH(resolution = resolution,
        #                        N = lengths[i],
        #                        alpha_par = alpha,
        #                        delta_par = delta,
        #                        p = 5)
        dat <- generateFunctionalCondVol_ff(
          resolution = resolution,
          N = lengths[i],
          alpha_par =
            paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\ARCH5send',
                   3,'p',p,'.rds'),
          delta_par = delta,
          burnin = 100,
          p = p)

        ######
        ## Estimation
        ######
        result_yw <- operatorARCH(X = dat$data, p = p )
        result_mp80 <- operatorARCHmp(X = dat$data, w = 0.2, p = p )

        alpha_true <- do.call(cbind, rep(list(dat$params$alpha),p))

        ## Get Errors
        errors_yw[sim] <-
          dot_integrate(dot_integrate_col( (result_yw$hat_alpha - alpha_true)^2 ))
        errors_yw1[sim] <-
          dot_integrate(dot_integrate_col(
            (result_yw$hat_alpha[,1:nrow(result_yw$hat_alpha)] -
               alpha_true[,1:nrow(result_yw$hat_alpha)] )^2 ))
        errors_yw2[sim] <-
          dot_integrate(dot_integrate_col(
            (result_yw$hat_alpha[,nrow(result_yw$hat_alpha)+1:nrow(result_yw$hat_alpha)] -
               alpha_true[,nrow(result_yw$hat_alpha)+1:nrow(result_yw$hat_alpha)] )^2 ))

        errors_mp[sim] <-
          dot_integrate(dot_integrate_col( (result_mp80$hat_alpha - alpha_true)^2 ))
        errors_mp1[sim] <-
          dot_integrate(dot_integrate_col(
            (result_mp80$hat_alpha[,1:nrow(result_mp80$hat_alpha)] -
               alpha_true[,1:nrow(result_mp80$hat_alpha)] )^2 ))
        errors_mp2[sim] <-
          dot_integrate(dot_integrate_col(
            (result_mp80$hat_alpha[,nrow(result_mp80$hat_alpha)+1:nrow(result_mp80$hat_alpha)] -
               alpha_true[,nrow(result_mp80$hat_alpha)+1:nrow(result_mp80$hat_alpha)] )^2 ))

      }
      cat(paste0('\n'))

      errors_list_yw[[i]] <- list(errors_yw,errors_yw1,errors_yw2)
      errors_list_mp[[i]] <- list(errors_mp,errors_mp1,errors_mp2)

      error_results_yw[i,] <- c(mean(errors_yw), var(errors_yw))
      error_results_yw1[i,] <- c(mean(errors_yw1), var(errors_yw1))
      error_results_yw2[i,] <- c(mean(errors_yw2), var(errors_yw2))
      error_results_mp[i,] <- c(mean(errors_mp), var(errors_mp))
      error_results_mp1[i,] <- c(mean(errors_mp1), var(errors_mp1))
      error_results_mp2[i,] <- c(mean(errors_mp2), var(errors_mp2))
    }

    list(error_results_yw, error_results_yw1, error_results_yw2, errors_list_yw,
         error_results_mp, error_results_mp1, error_results_mp2, errors_list_mp)
  }

# res_mp <- res
res <- op5_alphaComparision()
res[[1]]
