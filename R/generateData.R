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
                                      burnin=1000){

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
           sim_sigma2_matrix[,1] <- delta_par(times)
           sim_garch_matrix[,1] <- delta_par(times) * error_matrix[,1]

           for(j in 2:NBurn){
             # Variance
             for(i in 1:resolution){
               vector_op <- alpha_par(times[i], times) *
                 ((sim_garch_matrix[,(j-1)])^2)
               # sigma^2 = delta + \int alpha(t,s) X(s) ds, X=y[,j-1]^2
               sim_sigma2_matrix[i,j] <- delta_par(times[i]) +
                 .approx_int(as.matrix(vector_op))
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
           sim_sigma2_matrix[,1] <- delta_par(times)
           sim_garch_matrix[,1] <- delta_par(times) * error_matrix[,1]

           for(j in 2:NBurn){
             # Variance
             for(i in 1:resolution){
               vector_for_alpha_op <- alpha_par(times[i], times) *
                 (sim_garch_matrix[,(j-1)]^2)
               vector_for_beta_op <- beta_par(times[i], times) *
                 (sim_sigma2_matrix[,j-1])
               sim_sigma2_matrix[i,j] <- delta_par(times[i]) +
                 .approx_int(as.matrix(vector_for_alpha_op)) +
                 .approx_int(as.matrix(vector_for_beta_op))
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
                   delta_par=delta_par))
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
generateFunctionalCondVol_ff <- function(resolution, N, type = 'arch',
                                         alpha_par=NULL, delta_par=NULL, beta_par=NULL,
                                         burnin=1000){

  times <- seq(0,1,length.out=resolution)

  switch(type,
         arch = {
           NBurn <- burnin + N
           sim_sigma2_matrix <- array(NA, c(resolution, resolution, NBurn))
           sim_garch_matrix <- matrix(NA,nrow=resolution,ncol=NBurn)
           eps <- .generateFunctionalIID(resolution = resolution, N = NBurn)
           eps2 <- array(apply(eps, 2, tcrossprod),
                         dim=c(nrow(eps),nrow(eps),ncol(eps)))

           C_eps <- matrix(NA,resolution,resolution)
           for (i in 1:resolution){
             C_eps[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
           }
           # C_eps <- .estim_lag_cov(eps)

           if(is.null(alpha_par))
             alpha_par <- function(s,t,u,v){ s*(1-s)*t*(1-t)*u*(1-u)*v*(1-v) }
           if(is.null(delta_par))
             delta_par <- function(s,t){s*(1-s)*t*(1-t)}

           # Define variance and GARCH matrices
           sim_garch_matrix[,1] <- sapply(times,delta_par, t=times) %*% eps[,1]

           vector_op1 <- vector_op <-
             array(dim = c(resolution,resolution,resolution))
           for(i in 2:NBurn){
             prev_tensor <- sim_garch_matrix[,i-1] %*% t(sim_garch_matrix[,i-1])
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
             sim_sigma2_matrix[,,i] <-
               sapply(times,delta_par,t=times) + rowMeans(vector_op1, dims=2)

             #    See (66)
             Sigma2_svd <- svd(sim_sigma2_matrix[,,i])
             Sigma2_sqrt <- Sigma2_svd$u %*% diag(sqrt(Sigma2_svd$d)) %*% t(Sigma2_svd$v)
             # Formula (2)
             sim_garch_matrix[,i] <- Sigma2_sqrt %*% eps[,i]
           }



           ## PREVIOUS
           # for(j in 2:NBurn){
           #   # TODO:: Change to matrix mult
           #   for(j in 1:dims[1]){
           #     vector_op <- alpha[j, ] * sim_garch_matrix[,i-1] # X^2
           #     estimate[j,i] <- delta[j] + .approx_int(as.matrix(vector_op))
           #   }
           #   # Variance
           #   tensor_vals <- sim_garch_matrix[,(j-1)] %*% t(sim_garch_matrix[,(j-1)])
           #   for(i in 1:resolution){
           #     for(k in 1:resolution){
           #       vector_op <- sapply(times,alpha_par, u=times[i], v=times[k], s=times) %*% tensor_vals
           #       # sigma^2 = delta + \int alpha(u,v,t,s) X(u)X(v) du dv, X=y[,j-1]
           #       sim_sigma2_matrix[i,j] <- sim_sigma2_matrix[i,k,1] +
           #         .approx_int(as.matrix(vector_op))
           #     }
           #   }
           #
           #   # GARCH(1,0)/ARCH(1)
           #   sim_garch_matrix[,j] <- sqrt(sim_sigma2_matrix[,j]) *
           #     error_matrix[,j]
           # }
         },
         # garch = {
         #   NBurn <- burnin + N
         #   sim_sigma2_matrix <- sim_garch_matrix <-
         #     matrix(NA, resolution, NBurn)
         #   error_matrix <- .generateFunctionalIID(resolution = resolution, N = NBurn)
         #
         #   # Define the alpha, beta, and delta functions (if not given)
         #   if( is.null(alpha_par) )
         #     alpha_par <- function(t,s){ 12 * t * (1-t) * s * (1-s) }
         #   if( is.null(beta_par) )
         #     beta_par <- function(t,s){ 12 * t * (1-t) * s * (1-s) }
         #   if(is.null(delta_par))
         #     delta_par <- function(t){ 0.5* t * (1-t) }
         #
         #   # Define variance and GARCH matrices
         #   sim_sigma2_matrix[,1] <- delta_par(times)
         #   sim_garch_matrix[,1] <- delta_par(times) * error_matrix[,1]
         #
         #   for(j in 2:NBurn){
         #     # Variance
         #     for(i in 1:resolution){
         #       vector_for_alpha_op <- alpha_par(times[i], times) *
         #         (sim_garch_matrix[,(j-1)]^2)
         #       vector_for_beta_op <- beta_par(times[i], times) *
         #         (sim_sigma2_matrix[,j-1])
         #       sim_sigma2_matrix[i,j] <- delta_par(times[i]) +
         #         .approx_int(as.matrix(vector_for_alpha_op)) +
         #         .approx_int(as.matrix(vector_for_beta_op))
         #     }
         #
         #     # GARCH(1,1)
         #     sim_garch_matrix[,j] <- sqrt(sim_sigma2_matrix[,j]) *
         #       error_matrix[,j]
         #   }
         # },
         stop("Make type 'arch'"))

  list(garch_mat = sim_garch_matrix[, burnin + 1:N],
       sigma_mat = sim_sigma2_matrix[,, burnin + 1:N],
       error_matrix = eps[, burnin + 1:N],
       params=list(alpha_par=alpha_par,
                   beta_par=beta_par,
                   delta_par=delta_par))
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
    comat[i,] <- exp(-times[i]/2-times/2) * pmin(exp(times[i]), exp(times))
  }

  fiid <- MASS::mvrnorm(n = N,
                        mu = c(rep(0,resolution)),
                        Sigma = comat,
                        empirical = FALSE)

  t(fiid)
}
