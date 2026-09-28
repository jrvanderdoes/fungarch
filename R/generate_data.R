#' Generate functional GARCH data
#'
#' Simulate observations from the package's functional GARCH model using a
#' basis expansion and a specified covariance structure.
#'
#' @param N Number of observations to generate after burn-in.
#' @param resolution Number of evaluation points for each functional
#'   observation.
#' @param alpha_scale Coefficients used to construct the autoregressive
#'   operator.
#' @return A list containing simulated data, conditional covariance matrices,
#'   innovations, and the model parameters.
#' @export
#'
#' @examples
#' \dontrun{
#' generate_data(N = 10, resolution = 5)
#' }
generate_data <- function(N=750, resolution=20,
                           alpha_scale = c(0.7,0.7)){
  if (length(N) != 1L || !is.numeric(N) || !is.finite(N) ||
      N < 1 || N != floor(N)) {
    stop("N must be a single positive integer", call. = FALSE)
  }
  if (length(resolution) != 1L || !is.numeric(resolution) ||
      !is.finite(resolution) || resolution < 1 || resolution != floor(resolution)) {
    stop("resolution must be a single positive integer", call. = FALSE)
  }
  if (!is.numeric(alpha_scale) || length(alpha_scale) < 1L ||
      any(!is.finite(alpha_scale)) || length(alpha_scale) > resolution) {
    stop("alpha_scale must contain between one and resolution finite numeric values",
         call. = FALSE)
  }

  N <- as.integer(N)
  resolution <- as.integer(resolution)

  ## Defns
  C_eps <- .generate_ceps(resolution)
  eig_Ce <- eigen( C_eps )
  basis <- eig_Ce$vectors * sqrt(resolution)
  eigenvalues <- eig_Ce$values / resolution

  ##################
  ## SETUP ALPHA/DELTA
  ##################
  #   Alpha is based on tensor sum of `alpha_dimen` basis functions
  alpha <- matrix(0, nrow=resolution^2, ncol = resolution^2)
  for(k in 1:length(alpha_scale)){
    hkk <- as.vector(basis[,k] %o% basis[,k])
    alpha <- alpha + alpha_scale[k] * (hkk %o% hkk)
  }

  #   Delta is full C_eps
  delta <- C_eps


  ##################
  ## Generate Data
  ##################
  data <- generateFARCH_use(N = N, resolution = resolution,
                            alpha=alpha, delta=delta,
                            basis=basis,
                            burnin=100, p=1,
                            eigs=eigenvalues[1:length(alpha_scale)],
                            C_eps=C_eps, kl_errors=FALSE)

  data
}


#' Generate functional GARCH observations using a supplied model.
#'
#' This internal generator simulates functional observations after a burn-in
#' period using either a Karhunen--Loève innovation expansion or multivariate
#' normal innovations.
#'
#' @param N Number of observations to return after burn-in.
#' @param resolution Number of evaluation points per observation.
#' @param alpha Autoregressive operator. For `p = 1`, this is a square matrix;
#'   for `p > 1`, it is a list of `p` square matrices.
#' @param delta Baseline covariance contribution, as a square matrix.
#' @param basis Basis-function matrix used for KL innovations.
#' @param burnin Number of initial observations to discard.
#' @param p Number of autoregressive lags.
#' @param eigs Positive eigenvalues for KL innovations.
#' @param C_eps Innovation covariance matrix used when `kl_errors = FALSE`.
#' @param kl_errors If `TRUE`, generate innovations from `basis` and `eigs`;
#'   otherwise use `MASS::mvrnorm()` with `C_eps`.
#'
#' @return A list containing `data`, `sigma_mat`, `error_matrix`, and `params`.
#'
#' @noRd
generateFARCH_use <- function(N, resolution, alpha, delta, basis,
                              burnin=100, p=1, eigs=NULL,
                              C_eps=NULL, kl_errors=TRUE){
  if (length(N) != 1L || !is.numeric(N) || !is.finite(N) ||
      N < 1 || N != floor(N)) {
    stop("N must be a single positive integer", call. = FALSE)
  }
  if (length(resolution) != 1L || !is.numeric(resolution) ||
      !is.finite(resolution) || resolution < 1 || resolution != floor(resolution)) {
    stop("resolution must be a single positive integer", call. = FALSE)
  }
  if (length(burnin) != 1L || !is.numeric(burnin) || !is.finite(burnin) ||
      burnin < 0 || burnin != floor(burnin)) {
    stop("burnin must be a single nonnegative integer", call. = FALSE)
  }
  if (length(p) != 1L || !is.numeric(p) || !is.finite(p) ||
      p < 1 || p != floor(p)) {
    stop("p must be a single positive integer", call. = FALSE)
  }
  NBurn <- burnin + N
  if (p >= NBurn) {
    stop("p must be smaller than burnin + N", call. = FALSE)
  }
  if (length(kl_errors) != 1L || !is.logical(kl_errors) || is.na(kl_errors)) {
    stop("kl_errors must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.matrix(basis) || !is.numeric(basis) ||
      nrow(basis) != resolution || ncol(basis) < 1L ||
      any(!is.finite(basis))) {
    stop("basis must be a finite numeric matrix with resolution rows", call. = FALSE)
  }
  if (!is.matrix(delta) || !is.numeric(delta) ||
      nrow(delta) != resolution || ncol(delta) != resolution ||
      any(!is.finite(delta))) {
    stop("delta must be a finite resolution by resolution numeric matrix", call. = FALSE)
  }
  if (p == 1L) {
    if (!is.matrix(alpha) || !is.numeric(alpha) ||
        nrow(alpha) != resolution^2 || ncol(alpha) != resolution^2 ||
        any(!is.finite(alpha))) {
      stop("for p = 1, alpha must be a finite resolution^2 by resolution^2 matrix",
           call. = FALSE)
    }
  } else if (!is.list(alpha) || length(alpha) != p ||
             any(vapply(alpha, function(a) {
               !is.matrix(a) || !is.numeric(a) ||
                 nrow(a) != resolution^2 || ncol(a) != resolution^2 ||
                 any(!is.finite(a))
             }, logical(1)))) {
    stop("for p > 1, alpha must be a list of p finite resolution^2 by resolution^2 matrices",
         call. = FALSE)
  }
  if (!is.null(C_eps) &&
      (!is.matrix(C_eps) || !is.numeric(C_eps) ||
       nrow(C_eps) != resolution || ncol(C_eps) != resolution ||
       any(!is.finite(C_eps)))) {
    stop("C_eps must be a finite resolution by resolution numeric matrix", call. = FALSE)
  }
  if (kl_errors && (is.null(eigs) || length(eigs) == 0L)) {
    stop("eigs must be supplied when kl_errors = TRUE", call. = FALSE)
  }
  if (!kl_errors && is.null(C_eps)) {
    stop("C_eps must be supplied when kl_errors = FALSE", call. = FALSE)
  }
  if (!is.null(eigs) &&
      (!is.numeric(eigs) || any(!is.finite(eigs)) || any(eigs <= 0))) {
    stop("eigs must contain positive finite numeric values", call. = FALSE)
  }
  if (!is.null(eigs) && length(eigs) > ncol(basis)) {
    stop("eigs cannot contain more values than basis has columns", call. = FALSE)
  }

  #### Functions
  .get_kl_errors <- function(N, resolution, basis, eigs){

    .get_errors <- function(basis,eigs){
      ## NEXT
      xi <- sapply(eigs,
                   function(e) {
                     stats::rnorm(1, mean = 0, sd = sqrt(e))
                   }
      )

      Zeta <- tryCatch(basis[,1:length(eigs),drop=FALSE] *
                         xi[col(basis[,1:length(eigs),drop=FALSE])],
                       error = function(e) {
                         stop(call. = F, paste0(
                           "Check number of eigenvalues given. ",
                           "It does not match number of basis functions. ",
                           "Note, did you account for the constant function if ",
                           "it is in the basis?"
                         ))
                       }
      )

      rowSums(Zeta)
    }

    # OLD: data <- data.frame(matrix(NA, ncol = N, nrow = resolution))
    data <- matrix(NA_real_, ncol = N, nrow = resolution)
    for (j in 1:N) {
      data[, j] <- .get_errors(basis,eigs)
    }

    data
  }

  #### Code
  # OLD: times <- seq(0, 1, length.out=resolution)

  sim_sigma2_matrix2 <- matrix(0, ncol = NBurn, nrow = resolution^2)
  arch_data2 <- matrix(NA, nrow = resolution, ncol = NBurn)

  if(kl_errors){
    eps <- .get_kl_errors(N = NBurn, resolution = resolution,
                          basis=basis, eigs=eigs)
  }else{
    eps <- t(MASS::mvrnorm(n = NBurn,
                           mu = c(rep(0,resolution)),
                           Sigma = C_eps,
                           empirical = FALSE))
  }

  ## Define variance and GARCH matrices
  # Initialize
  for(i in 1:p){
    arch_data2[,i] <- operator(delta, eps[,i])
  }


  # Set matrices
  # OLD: Sigma2_list1 <- Sigma2_sqrt_list1 <-
  #   Sigma2_list2 <- Sigma2_sqrt_list2 <- list()
  # OLD: Sigma2_list2 <- list()
  # OLD: Sigma2_sqrt_list2 <- list()

  # Compute values
  for(i in (p+1):NBurn){
    tmp_alpha1 <- tmp_alpha2 <- rep(0, resolution^2)
    for(pp in 1:p){
      # Alpha
      prev_tensor <- tensor_vec(arch_data2[,i-pp], arch_data2[,i-pp])
      if(p==1){
        tmp_alpha2 <- tmp_alpha2 + operator(alpha, prev_tensor)
      }else{
        # tmp_alpha2 <- tmp_alpha2 + operator(alpha[[p]], prev_tensor)
        tmp_alpha2 <- tmp_alpha2 + operator(alpha[[pp]], prev_tensor)
      }
    }

    #  Sigma^2
    sim_sigma2_matrix2[,i] <- as.vector(delta) + tmp_alpha2

    # Compute the sqrt matrix
    tryCatch({
      Sigma2_svd2 <- svd(matrix(sim_sigma2_matrix2[,i],
                                ncol=resolution,nrow=resolution))
    }, error=function(e){
      stop('Error in sqrt matrix')
    })
    Sigma2_sqrt2 <- sqrt(resolution) *
      Sigma2_svd2$u %*% diag(sqrt(Sigma2_svd2$d)) %*% t(Sigma2_svd2$v)
    arch_data2[,i] <- operator(Sigma2_sqrt2, eps[,i])

    # OLD: Sigma2_sqrt_list2[[i]] <- Sigma2_sqrt2
  }

  # Return
  new_res1 <- list(
    data = arch_data2[, burnin + 1:N],
    sigma_mat = sim_sigma2_matrix2[, burnin + 1:N],
    error_matrix = eps[, burnin + 1:N],
    params=list(p=p, alpha = alpha, delta = delta)
  )
}
