
generateFunctionalARCH_both <- function(resolution, N,
                                         delta_par=NULL,
                                         burnin=100, p=1,
                                         K=2, scale=0.4){
  times <- seq(0, 1, length.out=resolution)

  NBurn <- burnin + N
  # sim_sigma2_matrix1 <-
  sim_sigma2_matrix2 <-
    matrix(0, ncol = NBurn, nrow = resolution^2)
  # arch_data <- arch_data1 <-
  arch_data2 <-
    matrix(NA, nrow = resolution, ncol = NBurn)

  eps <- .generateOU(resolution = resolution, N = NBurn)

  ## Define variance and GARCH matrices
  # Initialize
  for(i in 1:p){
    # arch_data1[,i] <- operator(sapply(times,delta_par, t=times), eps[,i])
    arch_data2[,i] <- operator(sapply(times,delta_par, t=times), eps[,i])
  }

  # alpha <- readRDS(alpha_par)

  # Set matrices
  Sigma2_list1 <- Sigma2_sqrt_list1 <-
    Sigma2_list2 <- Sigma2_sqrt_list2 <- list()

  # SETUP BASIS
  val_seq <- seq(0, 1, length.out=resolution)
  tmp <- fda::create.fourier.basis(nbasis = K, rangeval = c(0,1))
  basis_vals <- fda::eval.basis(seq(0,1,length.out=length(val_seq)), tmp)

  # ALPHA
  alpha <- matrix(0,nrow=resolution^2,ncol = resolution^2)
  for(k in 2:K){
    hkk <- tensor_vec(basis_vals[,k], basis_vals[,k])
    alpha <- alpha + tensor_vec(hkk, hkk)
  }
  alpha <- scale * alpha


  # Compute values
  for(i in (p+1):NBurn){
    tmp_alpha1 <- tmp_alpha2 <- rep(0, resolution^2)
    for(pp in 1:p){
      # Basis
      prev_tensor <- tensor_vec(arch_data1[,i-pp], arch_data1[,i-pp])
      # for(b in 2:K){
      #   mat_val <- basis_vals[,b] %o% basis_vals[,b]
      #   tmp <- mat_val * scale *
      #     as.numeric((basis_vals[,b] %*% arch_data1[,i-pp] / nrow(basis_vals))^2)
      #     # as.numeric((basis_vals[,b] %*% arch_data1[,i-pp])^2 / nrow(basis_vals)^2)
      #   tmp_alpha1 <- tmp_alpha1 + tmp
      # }

      # Alpha
      tmp_alpha2 <- tmp_alpha2 + operator(alpha, prev_tensor)
    }

    #  Sigma^2
    # sim_sigma2_matrix1[,i] <-
    #   as.vector(sapply(times,delta_par,t=times)) +
    #   tmp_alpha1
    sim_sigma2_matrix2[,i] <-
      as.vector(sapply(times,delta_par,t=times)) +
      tmp_alpha2

    # Attempt 2: Compute the sqrt matrix
    tryCatch({
      # Sigma2_svd1 <- svd(matrix(sim_sigma2_matrix1[,i],
      #                           ncol=resolution,nrow=resolution))
      Sigma2_svd2 <- svd(matrix(sim_sigma2_matrix2[,i],
                                ncol=resolution,nrow=resolution))
    }, error=function(e){
      cat(i)
      # View(sim_sigma2_matrix1)
      # View(arch_data1)
      stop('Error')
    })
    # Sigma2_sqrt1 <- sqrt(resolution) *
    #   Sigma2_svd1$u %*% diag(sqrt(Sigma2_svd1$d)) %*% t(Sigma2_svd1$v)
    Sigma2_sqrt2 <- sqrt(resolution) *
      Sigma2_svd2$u %*% diag(sqrt(Sigma2_svd2$d)) %*% t(Sigma2_svd2$v)
    # arch_data1[,i] <- operator(Sigma2_sqrt1, eps[,i])
    arch_data2[,i] <- operator(Sigma2_sqrt2, eps[,i])

    # Sigma2_sqrt_list1[[i]] <- Sigma2_sqrt1
    Sigma2_sqrt_list2[[i]] <- Sigma2_sqrt2
  }

  # Return
  new_res1 <- list(#data = arch_data1[, burnin + 1:N],
                   data = arch_data2[, burnin + 1:N],
                   # sigma_mat = sim_sigma2_matrix1[, burnin + 1:N],
                   sigma_mat = sim_sigma2_matrix2[, burnin + 1:N],
                   error_matrix = eps[, burnin + 1:N],
                   params=list(p=p,q=q,
                               # alpha_par=alpha_par,
                               alpha = alpha,
                               delta = sapply(times,delta_par,t=times))
  )
}

########################################################
# set.seed(123456)
dat_tmp <- generateFunctionalARCH_both(
  resolution = 20,
  N = 10000,
  delta_par = function(t,s){pmin(t,s)+0.01},
  burnin = 100,
  p = 1, K=2)
sqrt( sum(dat_tmp$params$alpha^2)/20^4 )
sqrt(.95^4*.4^2)
# should be scale

# # acf(colSums(dat_tmp$data^2))
# plot_fd(dat_tmp$data[,4800:5000])
alpha_true_tmp <- dat_tmp$params$alpha
plot_fd(alpha_true_tmp) # Max of 1.5ish

# Compare chat, dhat
result_mp95_tmp <- .est_opARCH_mp(X = dat_tmp$data1, w = 0.05, p = 1 )
result_yw95_tmp <- .operatorARCH(X = dat_tmp$data1, w = 0.05, p = 1 )
plot_fd(result_mp95_tmp$hat_alpha)
plot_fd(result_yw95_tmp$hat_alpha)


########################################################
set.seed(123456)
dat_base <- generateFunctionalARCH_both(
  resolution = 20,
  N = 750,
  delta_par = function(t,s){pmin(t,s)+0.01},
  burnin = 100,
  p = 1)
acf(colSums(dat_base$data^2))
acf(colSums(dat_base$data1^2))
plot_fd(dat_base$data)
sum(dat_base$data-dat_base$data1)

result_mp80_base <- .est_opARCH_mp(X = dat_base$data1, w = 0.2, p = 1 )
result_mp90_base <- .est_opARCH_mp(X = dat_base$data1, w = 0.1, p = 1 )
result_mp95_base <- .est_opARCH_mp(X = dat_base$data1, w = 0.05, p = 1 )
result_pw95_base <- .operatorARCH(X = dat_base$data1,p = 1)
result_pw999_base <- .operatorARCH(X = dat_base$data1,w = 0.001,p = 1)
alpha_true_base <- dat_base$params$alpha

dot_integrate(dot_integrate_col(
  (result_mp80_base$hat_alpha - alpha_true_base)^2 ))
dot_integrate(dot_integrate_col(
  (result_mp90_base$hat_alpha - alpha_true_base)^2 ))
dot_integrate(dot_integrate_col(
  (result_mp95_base$hat_alpha - alpha_true_base)^2 ))
dot_integrate(dot_integrate_col(
  (result_pw95_base$hat_alpha - alpha_true_base)^2 ))
dot_integrate(dot_integrate_col(
  (result_pw999_base$hat_alpha - alpha_true_base)^2 ))

plot_fd(alpha_true_base) # Max of 4ish
plot_fd(result_mp80_base$hat_alpha) # Way LARGE
plot_fd(result_mp90_base$hat_alpha) # Way LARGE
plot_fd(result_mp95_base$hat_alpha) # Way LARGE
plot_fd(result_pw95_base$hat_alpha) # Way small
plot_fd(result_pw999_base$hat_alpha) # 40 times low

###################################################
set.seed(123456)
dat <- generateFunctionalARCH_both(
  resolution = 40,
  N = 750,
  delta_par = function(t,s){pmin(t,s)+0.01},
  burnin = 100,
  p = 1)
acf(colSums(dat$data^2))
acf(colSums(dat$data1^2))
plot_fd(dat$data)
sum(dat$data-dat$data1)

result_mp80 <- .est_opARCH_mp(X = dat$data1, w = 0.2, p = 1 )
result_mp90 <- .est_opARCH_mp(X = dat$data1, w = 0.1, p = 1 )
result_mp95 <- .est_opARCH_mp(X = dat$data1, w = 0.05, p = 1 )
result_pw95 <- .operatorARCH(X = dat$data1,p = 1)
result_pw999 <- .operatorARCH(X = dat$data1,w = 0.001,p = 1)
alpha_true <- dat$params$alpha

dot_integrate(dot_integrate_col(
  (result_mp80$hat_alpha - alpha_true)^2 ))
dot_integrate(dot_integrate_col(
  (result_mp90$hat_alpha - alpha_true)^2 ))
dot_integrate(dot_integrate_col(
  (result_mp95$hat_alpha - alpha_true)^2 ))
dot_integrate(dot_integrate_col(
  (result_pw95$hat_alpha - alpha_true)^2 ))
dot_integrate(dot_integrate_col(
  (result_pw999$hat_alpha - alpha_true)^2 ))

plot_fd(alpha_true) # Max of 4ish
plot_fd(result_mp80$hat_alpha) # Way LARGE
plot_fd(result_mp90$hat_alpha) # Way LARGE
plot_fd(result_mp95$hat_alpha) # Way LARGE
plot_fd(result_pw95$hat_alpha) # 5x small
plot_fd(result_pw999$hat_alpha) # 5x small
#  max(alpha_true)/max(result_pw999$hat_alpha)

######################################

set.seed(1234567)
dat_long <- generateFunctionalARCH_both(
  resolution = 20,
  N = 1500,
  delta_par = function(t,s){pmin(t,s)+0.01},
  burnin = 100,
  p = 1)
acf(colSums(dat_long$data^2))
acf(colSums(dat_long$data1^2))
plot_fd(dat_long$data)
sum(dat_long$data-dat_long$data1)

result_mp80_long <- .est_opARCH_mp(X = dat_long$data1, w = 0.2, p = 1 )
result_mp90_long <- .est_opARCH_mp(X = dat_long$data1, w = 0.1, p = 1 )
result_mp95_long <- .est_opARCH_mp(X = dat_long$data1, w = 0.05, p = 1 )
result_pw95_long <- .operatorARCH(X = dat_long$data1,p = 1)
result_pw999_long <- .operatorARCH(X = dat_long$data1,w = 0.001,p = 1)
alpha_true_long <- dat_long$params$alpha

dot_integrate(dot_integrate_col(
  (result_mp80_long$hat_alpha - alpha_true_long)^2 ))
dot_integrate(dot_integrate_col(
  (result_mp90_long$hat_alpha - alpha_true_long)^2 ))
dot_integrate(dot_integrate_col(
  (result_mp95_long$hat_alpha - alpha_true_long)^2 ))
dot_integrate(dot_integrate_col(
  (result_pw95_long$hat_alpha - alpha_true_long)^2 ))
dot_integrate(dot_integrate_col(
  (result_pw999_long$hat_alpha - alpha_true_long)^2 ))

plot_fd(alpha_true_long) # Max of 4ish
plot_fd(result_mp80_long$hat_alpha) # Way LARGE
plot_fd(result_mp90_long$hat_alpha) # Way LARGE
plot_fd(result_mp95_long$hat_alpha) # Way LARGE
plot_fd(result_pw95_long$hat_alpha) #60x low
plot_fd(result_pw999_long$hat_alpha) # 60x times low


######################################

set.seed(1234567)
dat_low <- generateFunctionalARCH_both(
  resolution = 10,
  N = 750,
  delta_par = function(t,s){pmin(t,s)+0.01},
  burnin = 100,
  p = 1)
acf(colSums(dat_low$data^2))
acf(colSums(dat_low$data1^2))
plot_fd(dat_low$data)
sum(dat_low$data-dat_low$data1)

result_mp80_low <- .est_opARCH_mp(X = dat_low$data1, w = 0.2, p = 1 )
result_mp90_low <- .est_opARCH_mp(X = dat_low$data1, w = 0.1, p = 1 )
result_mp95_low <- .est_opARCH_mp(X = dat_low$data1, w = 0.05, p = 1 )
result_pw95_low <- .operatorARCH(X = dat_low$data1,p = 1)
result_pw999_low <- .operatorARCH(X = dat_low$data1,w = 0.001,p = 1)
alpha_true_low <- dat_low$params$alpha

dot_integrate(dot_integrate_col(
  (result_mp80_low$hat_alpha - alpha_true_low)^2 ))
dot_integrate(dot_integrate_col(
  (result_mp90_low$hat_alpha - alpha_true_low)^2 ))
dot_integrate(dot_integrate_col(
  (result_mp95_low$hat_alpha - alpha_true_low)^2 ))
dot_integrate(dot_integrate_col(
  (result_pw95_low$hat_alpha - alpha_true_low)^2 ))
dot_integrate(dot_integrate_col(
  (result_pw999_low$hat_alpha - alpha_true_low)^2 ))

plot_fd(alpha_true_low) # Max of 4ish
plot_fd(result_mp80_low$hat_alpha) # Wrong shape
plot_fd(result_mp90_low$hat_alpha) # Wrong shape
plot_fd(result_mp95_low$hat_alpha) # Wrong shape
plot_fd(result_pw95_low$hat_alpha) # Wrong shape
plot_fd(result_pw999_low$hat_alpha) # Wrong shape, 230,000 too low


######################################################################
dat_K6 <- generateFunctionalARCH_both(
  resolution = 20,
  N = 750,
  delta_par = function(t,s){pmin(t,s)+0.01},
  burnin = 100,
  p = 1,
  K = 6)
acf(colSums(dat_K6$data^2))
acf(colSums(dat_K6$data1^2))
plot_fd(dat_K6$data)
sum(dat_K6$data-dat_K6$data1)

result_mp80_K6 <- .est_opARCH_mp(X = dat_K6$data1, w = 0.2, p = 1 )
result_mp90_K6 <- .est_opARCH_mp(X = dat_K6$data1, w = 0.1, p = 1 )
result_mp95_K6 <- .est_opARCH_mp(X = dat_K6$data1, w = 0.05, p = 1 )
result_pw95_K6 <- .operatorARCH(X = dat_K6$data1,p = 1)
result_pw999_K6 <- .operatorARCH(X = dat_K6$data1,w = 0.001,p = 1)
alpha_true_K6 <- dat_K6$params$alpha

dot_integrate(dot_integrate_col(
  (result_mp80_K6$hat_alpha - alpha_true_K6)^2 ))
dot_integrate(dot_integrate_col(
  (result_mp90_K6$hat_alpha - alpha_true_K6)^2 ))
dot_integrate(dot_integrate_col(
  (result_mp95_K6$hat_alpha - alpha_true_K6)^2 ))
dot_integrate(dot_integrate_col(
  (result_pw95_K6$hat_alpha - alpha_true_K6)^2 ))
dot_integrate(dot_integrate_col(
  (result_pw999_K6$hat_alpha - alpha_true_K6)^2 ))

plot_fd(alpha_true_K6) # Max of 11ish
plot_fd(result_mp80_K6$hat_alpha) # Wrong shape
plot_fd(result_mp90_K6$hat_alpha) # Wrong shape
plot_fd(result_mp95_K6$hat_alpha) # Wrong shape
plot_fd(result_pw95_K6$hat_alpha) # Wrong shape
plot_fd(result_pw999_K6$hat_alpha) # Wrong shape, 480 too low


# K=1
# reduce alpha
