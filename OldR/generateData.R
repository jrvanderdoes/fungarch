#' Data Generating Process - Functional ARCH/GARCH Process
#'
#' @description Generates functional curve data following
#'  the functional ARCH(1) or GARCH(1,1) process.
#'
#' @param resolution The number of grid point in each curve observation.
#' @param N The sample size.
#' @param type A string to switch data generating process between "arch" -
#'  functional ARCH and "garch" - functional GARCH.
#' @param alpha_par The ARCH kernel coefficient function in the conditional
#'  volatility equation. If it is missing, "\eqn{12 * t * (1-t) * s * (1-s)}"
#'  is used to generate FGARCH, and "\eqn{16 * t * (1-t) * s * (1-s)}" is used
#'  to generate FARCH, for \eqn{t\in[0,1]} and \eqn{s\in[0,1]}.
#' @param beta_par The GARCH kernel coefficient function in the conditional
#'  volatility equation. If it is missing, "\eqn{12 * t * (1-t) * s * (1-s)}" is
#'  used to generate FGARCH, for \eqn{t\in[0,1]} and \eqn{s\in[0,1]}.
#'
#' @return List of generated processes:
#' @return garch_mat: FARCH/GARCH sequences, where the finite realization of
#'  curves are stored in columns;
#' @return sigma_mat: Conditional volatility sequences,  where the finite
#'  realization of curves are stored in columns.
#'
#' @export
#'
#' @details
#' If \eqn{x_i(t)} follows an FARCH(1) process,\cr
#' \eqn{x_i(t)=\sigma_i(t)\varepsilon_i(t)}, \eqn{t \in [0,1]},\cr
#' \eqn{\sigma_i^2(t)=\omega(t)+\int \alpha(t,s) x^2_{i-1}(s)ds}.\cr
#'
#' If \eqn{x_i(t)} follows an FGARCH(1,1) process,\cr
#' \eqn{x_i(t)=\sigma_i(t)\varepsilon_i(t)}, \eqn{t \in [0,1]},\cr
#' \eqn{\sigma_i^2(t)=\omega(t)+\int \alpha(t,s) x^2_{i-1}(s)ds+\int \beta(t,s) \sigma^2_{i-1}(s)ds},\cr
#' where the innovation \eqn{\varepsilon_i(t)} follows an Ornstein–Uhlenbeck process \code{\link{.generateFunctionalIID}}, and the constant coefficient \eqn{\omega(t)=0.1t(1-t)}.
#'
#' @seealso \code{\link{.generateFunctionalIID}}
#' @examples
#' # generate discrete evaluations of 100 FARCH curves that each curve is
#'  realized on 50 grid points.
#' data <- generateFunctionalCondVol(50,100,'garch')
#' plot_fd(data$garch_mat)
#' data <- generateFunctionalCondVol(50,100,'arch')
#' plot_fd(data$garch_mat)
#'
#' @references
#' Hormann, S., Horvath, L., Reeder, R. (2013). A functional version of the ARCH model. Econometric Theory. 29(2), 267-288. <doi:10.1017/S0266466612000345>.\cr
#' Aue, A., Horvath, L., F. Pellatt, D. (2017). Functional generalized autoregressive conditional heteroskedasticity. Journal of Time Series Analysis. 38(1), 3-21. <doi:10.1111/jtsa.12192>.\cr
#'
generateFunctionalCondVol <- function(resolution, N, type = 'garch',
                                      alpha_par=NULL, beta_par=NULL, delta_par=NULL,
                                      burnin=1000, p = 1, q=1){

  times <- seq(0,1,length.out=resolution)
  #times <- 1:resolution/resolution

  switch(type,
         arch = {
           NBurn <- burnin + N
           sim_sigma2_matrix <- sim_garch_matrix <-
             matrix(NA, resolution, NBurn)
           error_matrix <- .generateFunctionalIID(resolution = resolution, N = NBurn)
           #error_matrix <- matrix(1,nrow = resolution, ncol = NBurn)

           if(is.null(alpha_par))
             alpha_par <- function(t,s){ 16 * t * (1-t) * s * (1-s) }
           if(is.null(delta_par))
             delta_par <- function(t){ 0.1* t * (1-t) }

           # Define variance and GARCH matrices
           sim_sigma2_matrix[,1:p] <- delta_par(times)
           sim_garch_matrix[,1:p] <- delta_par(times) * error_matrix[,1]

           for(j in (p+1):NBurn){
             # Variance
             for(i in 1:resolution){
               # sigma^2 = delta + \int alpha(t,s) X(s) ds, X=y[,j-1]^2
               alpha_sum <- 0
               for(pp in 1:p){
                 vector_op <- alpha_par(times[i], times) *
                   ((sim_garch_matrix[,(j-pp)])^2)
                 alpha_sum <- alpha_sum +
                   dot_integrate(as.matrix(vector_op))
               }
               sim_sigma2_matrix[i,j] <- delta_par(times[i]) +
                 alpha_sum
             }

             # GARCH(1,0)/ARCH(1)
             sim_garch_matrix[,j] <- sqrt(sim_sigma2_matrix[,j]) *
               error_matrix[,j]
           }
         },
         garch = {
           NBurn <- burnin + N
           sim_sigma2_matrix <- sim_garch_matrix <-
             matrix(NA, resolution, NBurn)
           error_matrix <- .generateFunctionalIID(resolution = resolution, N = NBurn)

           # Define the alpha, beta, and delta functions (if not given)
           if( is.null(alpha_par) )
             alpha_par <- function(t,s){ 12 * t * (1-t) * s * (1-s) }
           if( is.null(beta_par) )
             beta_par <- function(t,s){ 12 * t * (1-t) * s * (1-s) }
           if(is.null(delta_par))
             delta_par <- function(t){ 0.5* t * (1-t) }

           # Define variance and GARCH matrices
           sim_sigma2_matrix[,1:q] <- delta_par(times)
           sim_garch_matrix[,1:p] <- delta_par(times) * error_matrix[,1]

           for(j in 2:NBurn){
             # Variance
             for(i in 1:resolution){
               alpha_sum <- 0
               for(qq in 1:p){
                 vector_for_alpha_op <- alpha_par(times[i], times) *
                   (sim_garch_matrix[,(j-qq)]^2)
                 alpha_sum <- alpha_sum +
                   dot_integrate(as.matrix(vector_for_alpha_op))
               }
               beta_sum <- 0
               for(qq in 1:q){
                 vector_for_beta_op <- beta_par(times[i], times) *
                   (sim_sigma2_matrix[,j-qq])
                 beta_sum <- beta_sum +
                   dot_integrate(as.matrix(vector_for_beta_op))
               }

               sim_sigma2_matrix[i,j] <- delta_par(times[i]) +
                 alpha_sum + beta_sum

             }

             # GARCH(1,1)
             sim_garch_matrix[,j] <- sqrt(sim_sigma2_matrix[,j]) *
               error_matrix[,j]
           }
         },
         stop("Make type 'arch' or 'garch'"))

  list(garch_mat = sim_garch_matrix[, burnin + 1:N],
       sigma_mat = sim_sigma2_matrix[, burnin + 1:N],
       error_matrix = error_matrix[, burnin + 1:N],
       params=list(alpha_par=alpha_par,
                   beta_par=beta_par,
                   delta_par=delta_par,
                   p=p, q=q))
}

#' Data Generating Process - Functional ARCH/GARCH Process
#'
#' @description Generates functional curve data following
#'  the functional ARCH(1) or GARCH(1,1) process.
#'
#' @param resolution The number of grid point in each curve observation.
#' @param N The sample size.
#' @param type A string to switch data generating process between "arch" -
#'  functional ARCH and "garch" - functional GARCH.
#' @param alpha_par The ARCH kernel coefficient function in the conditional
#'  volatility equation. If it is missing, "\eqn{12 * t * (1-t) * s * (1-s)}"
#'  is used to generate FGARCH, and "\eqn{16 * t * (1-t) * s * (1-s)}" is used
#'  to generate FARCH, for \eqn{t\in[0,1]} and \eqn{s\in[0,1]}.
#' @param beta_par The GARCH kernel coefficient function in the conditional
#'  volatility equation. If it is missing, "\eqn{12 * t * (1-t) * s * (1-s)}" is
#'  used to generate FGARCH, for \eqn{t\in[0,1]} and \eqn{s\in[0,1]}.
#'
#' @return List of generated processes:
#' @return garch_mat: FARCH/GARCH sequences, where the finite realization of
#'  curves are stored in columns;
#' @return sigma_mat: Conditional volatility sequences,  where the finite
#'  realization of curves are stored in columns.
#'
#' @export
#'
#' @details
#' If \eqn{x_i(t)} follows an FARCH(1) process,\cr
#' \eqn{x_i(t)=\sigma_i(t)\varepsilon_i(t)}, \eqn{t \in [0,1]},\cr
#' \eqn{\sigma_i^2(t)=\omega(t)+\int \alpha(t,s) x^2_{i-1}(s)ds}.\cr
#'
#' If \eqn{x_i(t)} follows an FGARCH(1,1) process,\cr
#' \eqn{x_i(t)=\sigma_i(t)\varepsilon_i(t)}, \eqn{t \in [0,1]},\cr
#' \eqn{\sigma_i^2(t)=\omega(t)+\int \alpha(t,s) x^2_{i-1}(s)ds+\int \beta(t,s) \sigma^2_{i-1}(s)ds},\cr
#' where the innovation \eqn{\varepsilon_i(t)} follows an Ornstein–Uhlenbeck process \code{\link{.generateFunctionalIID}}, and the constant coefficient \eqn{\omega(t)=0.1t(1-t)}.
#'
#' @seealso \code{\link{.generateFunctionalIID}}
#' @examples
#' data <- generateFunctionalCondVol_ff(24,100,'arch')
#'
#' @references
#'
generateFunctionalCondVol_ff <- function(resolution, N, type = 'arch',
                                         alpha_par=NULL, delta_par=NULL,
                                         beta_par=NULL, burnin=100,
                                         p=1, q=1){
  times <- seq(0,1,length.out=resolution)

  switch(type,
         arch = {
           # NBurn <- burnin + N
           # sim_sigma2_matrix <- array(NA, c(resolution, resolution, NBurn))
           # sim_garch_matrix <- matrix(NA,nrow=resolution,ncol=NBurn)
           # eps <- .generateFunctionalIID(resolution = resolution, N = NBurn)
           # eps2 <- array(apply(eps, 2, tcrossprod),
           #               dim=c(nrow(eps),nrow(eps),ncol(eps)))
           #
           # if(is.null(alpha_par))
           #   alpha_par <- function(s,t,u,v){ s*(1-s)*t*(1-t)*u*(1-u)*v*(1-v) }
           # if(is.null(delta_par))
           #   delta_par <- function(s,t){s*(1-s)*t*(1-t)}
           #
           # # Define variance and GARCH matrices
           # sim_garch_matrix[,1:p] <-
           #   sapply(times,delta_par, t=times) %*% eps[,1] / resolution
           #
           # vector_op1 <- vector_op <-
           #     array(dim = c(resolution,resolution,resolution))
           # # Try For speed
           # vector_list <- list()
           # for(j in 1:resolution){
           #   vector_list[[j]] <- list()
           #   # Over all t
           #   for(k in 1:resolution){
           #     # Makes first index of list s
           #     #  Makes second index of list t
           #     #  Makes rows of each u
           #     #  Makes cols of each v
           #     vector_list[[j]][[k]] <-
           #       as.matrix(
           #         sapply(times, alpha_par,
           #                s=times[j],t=times[k],u=times)
           #       )
           #   }
           # }
           #
           # # idx <- 0
           # # alpha <- matrix(nrow=resolution^2,ncol=resolution^2)
           # # for(j in 1:resolution){
           # #   for(k in 1:resolution){
           # #     idx <- idx + 1
           # #     alpha[,idx] <- as.vector(vector_list[[j]][[k]])
           # #   }
           # # }
           #
           # Sigma2_list <- Sigma2_sqrt_list <- list()
           #
           # for(i in (p+1):NBurn){
           #   tmp_alpha <- array(0,dim = c(resolution,resolution))
           #   # use_alpha <- rep(0, resolution^2)
           #   for(pp in 1:p){
           #     prev_tensor <- sim_garch_matrix[,i-pp] %*%
           #       t(sim_garch_matrix[,i-pp])
           #
           #     # use_alpha <- use_alpha + alpha %*% as.vector(prev_tensor)
           #     for(j in 1:resolution){
           #       # Over all t
           #       for(k in 1:resolution){
           #         # vector_op[,,k] <-
           #         #   as.matrix(
           #         #     sapply(times, alpha_par,
           #         #            s=times[j],t=times[k],u=times)
           #         #   ) %*% prev_tensor
           #         vector_op[,,k] <- vector_list[[j]][[k]] %*% prev_tensor / resolution
           #       }
           #       # TODO:: Convert to trap int
           #       vector_op1[,,j] <- rowMeans(vector_op, dims=2)
           #       # vector_op1_tmp[,,j] <- (vector_op[,,1] + vector_op[,,dim(vector_op)[3]] +
           #       #           rowSums(2*vector_op[,,2:(dim(vector_op)[3]-1)],dims=2) ) /
           #       #   (2*dim(vector_op)[3])
           #     }
           #     # Formula (2)
           #     #  delta + sum alpha(X2) + sum beta(Sigma2)
           #     # TODO:: Convert to trap int
           #     tmp_alpha <- tmp_alpha + rowMeans(vector_op1, dims=2)
           #   }
           #
           #
           #   sim_sigma2_matrix[,,i] <-
           #     sapply(times,delta_par,t=times) +
           #     tmp_alpha #rowMeans(vector_op1, dims=2)
           #
           #   #    See (66)
           #   Sigma2_svd <- svd(sim_sigma2_matrix[,,i])
           #   Sigma2_sqrt <- sqrt(resolution) * Sigma2_svd$u %*%
           #     diag(sqrt(Sigma2_svd$d)) %*%
           #     t(Sigma2_svd$v)
           #   # Formula (2)
           #   sim_garch_matrix[,i] <- Sigma2_sqrt %*% eps[,i] / resolution
           #
           #   Sigma2_sqrt_list[[i]] <- Sigma2_sqrt
           #   Sigma2_list[[i]] <- sim_sigma2_matrix[,,i]
           # }

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

           # # # Attempt 1: Precompute Alpha parameters to use as needed
           # # vector_op1 <- vector_op <-
           # #   array(dim = c(resolution,resolution,resolution))
           # # Try For speed
           # vector_list <- list()
           # for(j in 1:resolution){
           #   vector_list[[j]] <- list()
           #   # Over all t
           #   for(k in 1:resolution){
           #     # Makes first index of list s
           #     #  Makes second index of list t
           #     #  Makes rows of each u
           #     #  Makes cols of each v
           #     vector_list[[j]][[k]] <-
           #       as.matrix(
           #         sapply(times, alpha_par,
           #                s=times[j],t=times[k],u=times)
           #       )
           #   }
           # }
           #
           # # Attempt 2: Compute Alpha matrix
           # idx <- 0
           # alpha <- matrix(nrow=resolution^2,ncol=resolution^2)
           # for(j in 1:resolution){
           #   for(k in 1:resolution){
           #     idx <- idx + 1
           #     alpha[,idx] <- as.vector(t(vector_list[[j]][[k]]))
           #   }
           # }
           alpha <- readRDS(alpha_par)
             # paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\',
             #        'ARCH5send',3,'.rds'))

           # Set matrices
           # Sigma2_list <- Sigma2_sqrt_list <- list()
           Sigma2_list1 <- Sigma2_sqrt_list1 <- list()

           # Compute values
           for(i in (p+1):NBurn){
             # tmp_alpha <- array(0,dim = c(resolution,resolution))
             tmp_alpha1 <- rep(0, resolution^2)
             for(pp in 1:p){
               # array(apply(arch_data1[,i-pp,drop=FALSE], 2, tcrossprod),
               #       dim=c(nrow(arch_data1),nrow(arch_data1)))
               prev_tensor <-
                 as.vector( arch_data1[,i-pp] %*% t(arch_data1[,i-pp]) ) #/ resolution

               # Attempt 2:
               tmp_alpha1 <- tmp_alpha1 + alpha %*% prev_tensor / resolution^2

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
               cat(i)
               View(sim_sigma2_matrix1)
               View(arch_data1)
               stop('Error')
             })
             Sigma2_sqrt1 <- sqrt(resolution) * Sigma2_svd1$u %*%
               diag(sqrt(Sigma2_svd1$d)) %*% t(Sigma2_svd1$v)
             arch_data1[,i] <- Sigma2_sqrt1 %*% eps[,i] / resolution

             Sigma2_sqrt_list1[[i]] <- Sigma2_sqrt1
           }

           # Redefine later
           return(
             list(data = arch_data1[, burnin + 1:N],
                  sigma_mat = sim_sigma2_matrix1[, burnin + 1:N],
                  error_matrix = eps[, burnin + 1:N],
                  params=list(p=p,q=q,
                              alpha_par=alpha_par,
                              alpha = alpha,
                              delta_par=delta_par,
                              delta = as.vector(sapply(times,delta_par,t=times)),
                              Sigma2_list= Sigma2_list1[burnin + 1:N],
                              Sigma2_sqrt_list= Sigma2_sqrt_list1[burnin + 1:N]))
           )
         },
         garch = {
           NBurn <- burnin + N
           sim_sigma2_matrix <- array(0, c(resolution, resolution, NBurn))
           sim_garch_matrix <- matrix(NA,nrow=resolution,ncol=NBurn)
           eps <- .generateFunctionalIID(resolution = resolution, N = NBurn)
           eps2 <- array(apply(eps, 2, tcrossprod),
                         dim=c(nrow(eps),nrow(eps),ncol(eps)))

           if(is.null(alpha_par))
             alpha_par <- function(s,t,u,v){ s*(1-s)*t*(1-t)*u*(1-u)*v*(1-v) }
           if(is.null(beta_par))
             beta_par <- function(s,t,u,v){ s*(1-s)*t*(1-t)*u*(1-u)*v*(1-v) }
           if(is.null(delta_par))
             delta_par <- function(s,t){s*(1-s)*t*(1-t)}

           # Define variance and GARCH matrices
           sim_garch_matrix[,1] <-
             sapply(times,delta_par, t=times) %*% eps[,1] / resolution

           # Try For speed
           vector_list_alpha <- vector_list_beta <- list()
           for(j in 1:resolution){
             vector_list_alpha[[j]] <- list()
             vector_list_beta[[j]] <- list()
             # Over all t
             for(k in 1:resolution){
               # Makes first index of list s
               #  Makes second index of list t
               #  Makes rows of each u
               #  Makes cols of each v
               vector_list_alpha[[j]][[k]] <-
                 as.matrix(
                   sapply(times, alpha_par,
                          s=times[j],t=times[k],u=times)
                 )
               vector_list_beta[[j]][[k]] <-
                 as.matrix(
                   sapply(times, beta_par,
                          s=times[j],t=times[k],u=times)
                 )
             }
           }

           Sigma2_list <- Sigma2_sqrt_list <- list()
           vector_opA1 <- vector_opA <-
             vector_opB1 <- vector_opB <-
             array(dim = c(resolution,resolution,resolution))
           tmp_alpha <- tmp_beta <-
             array(0,dim = c(resolution,resolution))

           # SIMULATE DATA
           for(i in (max(p,q)+1):NBurn){
             # Check previous
             for(pq in 1:max(p,q)){
               prev_tensor <- sim_garch_matrix[,i-pq] %*% t(sim_garch_matrix[,i-pq])

               for(j in 1:resolution){
                 # Over all t
                 for(k in 1:resolution){
                   vector_opA[,,k] <-
                     vector_list_alpha[[j]][[k]] %*% (prev_tensor) / resolution
                   vector_opB[,,k] <-
                     vector_list_beta[[j]][[k]] %*% (sim_sigma2_matrix[,,i-pq]) / resolution
                 }
                 # TODO:: Convert to trap int
                 vector_opA1[,,j] <- rowMeans(vector_opA, dims=2)
                 vector_opB1[,,j] <- rowMeans(vector_opB, dims=2)
               }
               # Formula (2)
               #  delta + sum alpha(X2) + sum beta(Sigma2)
               # TODO:: Convert to trap int
               tmp_alpha <- rowMeans(vector_opA1, dims=2)
               tmp_beta <- rowMeans(vector_opB1, dims=2)
               if(pq<=p){
                 tmp_alpha <- tmp_alpha + rowMeans(vector_opA1, dims=2)
               }
               if(pq<=q){
                 tmp_beta <- tmp_beta + rowMeans(vector_opB1, dims=2)
               }
             }

             sim_sigma2_matrix[,,i] <-
               sapply(times,delta_par,t=times) +
               tmp_alpha + # rowMeans(vector_opA1, dims=2) +
               tmp_beta # rowMeans(vector_opB1, dims=2)

             #    See (66)
             Sigma2_svd <- svd(sim_sigma2_matrix[,,i])
             Sigma2_sqrt <- sqrt(resolution) * Sigma2_svd$u %*%
               diag(sqrt(Sigma2_svd$d)) %*%
               t(Sigma2_svd$v)
             # Formula (2)
             sim_garch_matrix[,i] <- Sigma2_sqrt %*% eps[,i] / resolution

             Sigma2_sqrt_list[[i]] <- Sigma2_sqrt
             Sigma2_list[[i]] <- sim_sigma2_matrix[,,i]
           }
           vector_list <- list(vector_opA1,vector_opB1)
         },
         stop("Make type 'arch' or 'garch'"))

  list(garch_mat = sim_garch_matrix[, burnin + 1:N],
       sigma_mat = sim_sigma2_matrix[,, burnin + 1:N],
       error_matrix = eps[, burnin + 1:N],
       params=list(p=p,q=q,
                   alpha_par=alpha_par,
                   alpha = alpha,
                   beta_par=beta_par,
                   delta_par=delta_par,
                   delta = as.vector(sapply(times,delta_par,t=times)),
                   Sigma2_list=
                     Sigma2_list[burnin + 1:N],
                   Sigma2_sqrt_list=
                     Sigma2_sqrt_list[burnin + 1:N],
                   vector_list=vector_list))
}


#' Data Generating Process - Independent Process
#'
#' @description Generates iid functional curve data following the
#' Ornstein–Uhlenbeck (OU) process.
#'
#' @param resolution The number of grid point in each curve observation.
#' @param N The sample size.
#'
#' @return A (resolution) x (number of observations) matrix for iid sequences,
#'  where the finite realization of curves are stored in columns.
#'
#' @details
#' \eqn{x_i(t)=e^{-t/2}W_i(e^t)}, \eqn{t \in [0,1]},\cr
#' where \eqn{W_i(t)} is a standard Brownian Motion.
#'
#' @seealso \code{\link{generateFunctionalCondVol}}
#' @examples
#' # generate discrete evaluations of 100 iid curves
#' #   each curve is realized on 50 grid points.
#' yd_iid <- .generateFunctionalIID(50, 100)
.generateFunctionalIID <- function(resolution, N){

  times <- 1:resolution/resolution

  # Covariance structure (OU process Cov)
  comat <- matrix(NA,resolution,resolution)
  for (i in 1:resolution){
    #comat[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
    comat[i,] <- exp(-abs(times[i]-times)/2)
  }

  fiid <- MASS::mvrnorm(n = N,
                        mu = c(rep(0,resolution)),
                        Sigma = comat,
                        empirical = FALSE)

  t(fiid)
}
