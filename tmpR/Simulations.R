#############################################
###         ARCH(1)
#############################################
alpha = function(s,t,u,v){
  val_seq <- seq(0,1,length.out=20)
  K <- 30 #30 #20 #15 #10 #5

  # tmp <- fda::create.fourier.basis(nbasis = K,rangeval = c(0,1))
  # basis_vals <- fda::eval.basis(seq(0,1,length.out=length(val_seq)), tmp)
  #
  # alp <- matrix(0,nrow=length(val_seq)^2,ncol = length(val_seq)^2)
  # for(k in 1:K){
  #   hkk <- basis_vals[,k] %*% t(basis_vals[,k]) /
  #     nrow(basis_vals)
  #   # tmp <- max(hkk)
  #   alp <- alp + 50*kronecker(hkk, hkk) / nrow(hkk)
  # }
  # alp <- (3.5*K)*(alp / K)*20
  # saveRDS(alp,
  #         paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'.rds'))
  alp <- readRDS(
    paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'.rds'))
  alp[(which(s==val_seq)-1)*length(val_seq)+which(t==val_seq),
      (which(u==val_seq)-1)*length(val_seq)+which(v==val_seq)]
}
delta = function(s,t){pmin(t,s)*20 }

dat <- generateFunctionalCondVol_ff(resolution = 20,
                                    N =750,type = 'arch',
                                    alpha_par = alpha,
                                    delta_par = delta)

acf(colSums(dat$garch_mat^2))
plot_fd(dat$garch_mat)

#############################################
###         GARCH(1,1)
############################################
alpha = function(s,t,u,v){
  val_seq <- seq(0,1,length.out=20)
  K <- 5 #30 #20 #15 #10 #5

  # tmp <- fda::create.fourier.basis(nbasis = K,rangeval = c(0,1))
  # basis_vals <- fda::eval.basis(seq(0,1,length.out=length(val_seq)), tmp)
  #
  # alp <- matrix(0,nrow=length(val_seq)^2,ncol = length(val_seq)^2)
  # for(k in 1:K){
  #   hkk <- basis_vals[,k] %*% t(basis_vals[,k]) /
  #     nrow(basis_vals)
  #   # tmp <- max(hkk)
  #   alp <- alp + 50*kronecker(hkk, hkk) / nrow(hkk)
  # }
  # alp <- 250 * alp / K
  # saveRDS(alp,
  #         paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\GARCHalp',K,'.rds'))
  alp <- readRDS(
    paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\GARCHalp',K,'.rds'))
  alp[(which(s==val_seq)-1)*length(val_seq)+which(t==val_seq),
      (which(u==val_seq)-1)*length(val_seq)+which(v==val_seq)]
}
beta = function(s,t,u,v){
  val_seq <- seq(0,1,length.out=20)
  K <- 5 #30 #20 #15 #10 #5

  # tmp <- fda::create.fourier.basis(nbasis = K,rangeval = c(0,1))
  # basis_vals <- fda::eval.basis(seq(0,1,length.out=length(val_seq)), tmp)
  #
  # alp <- matrix(0,nrow=length(val_seq)^2,ncol = length(val_seq)^2)
  # for(k in 1:K){
  #   hkk <- basis_vals[,k] %*% t(basis_vals[,k]) /
  #     nrow(basis_vals)
  #   # tmp <- max(hkk)
  #   alp <- alp + 50*kronecker(hkk, hkk) / nrow(hkk)
  # }
  # alp <- 250 * alp / K
  # saveRDS(alp,
  #         paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\GARCHbet',K,'.rds'))
  alp <- readRDS(
    paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\GARCHbet',K,'.rds'))
  alp[(which(s==val_seq)-1)*length(val_seq)+which(t==val_seq),
      (which(u==val_seq)-1)*length(val_seq)+which(v==val_seq)]
}
delta = function(s,t){ 50*pmin(t,s)  }

dat <- generateFunctionalCondVol_ff(resolution = 20,
                                    N = 750,'garch',
                                    alpha_par = alpha,
                                    beta_par = beta,
                                    delta_par = delta)
acf(colSums(dat$garch_mat^2))
plot_fd(dat$garch_mat)
###################################################
###         Pointwise - ARCH(1)
###################################################
library(devtools);load_all()
nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'YW'=NA,
             'MP90'=NA,
             'MP95'=NA,
             'MP99'=NA,
             'PW'=NA,
             'MLE'=NA,
             'Mean'=NA,'Prev'=NA,'Zero'=NA)

cat('Fit (',nSims,'): ')
set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitComparePWARCH(N_fit = 500, N_pred = 250,
                             alpha = function(t,s){ 10 * t * (1-t) * s * (1-s) },
                             delta = function(t){ rep(0.5,length(t)) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}
tmp <- round(colMeans(results_oneStep,na.rm = T),3)
round((tmp / min(tmp) - 1) * 100,3)

# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})
# apply(results_oneStep,MARGIN = 2,
#       function(X){
#         x = X[order(-X)]
#         round(mean(x[1:(0.95*length(x))]),3)
#       })
saveRDS(results_oneStep,paste0(getwd(),'/paperData/','pw.rds'))

###################################################
###         Pointwise - ARCH(5)
###################################################

###################################################
###         Pointwise - GARCH(1,1)
###################################################
library(devtools);load_all()
nSims <- 5#100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'YW'=0,
             'MP90'=0, 'MP95'=0, 'MP99'=0,
             'PW'=0,
             'MLE'=0,
             'Mean'=0, 'Prev'=0, 'Zero'=0)
cat('Fit (',nSims,'): ')
set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitComparePWGARCH(N_fit = 500, N_pred = 250,
                             alpha = function(t,s){ 10 * t * (1-t) * s * (1-s) },
                             beta = function(t,s){ 10 * t * (1-t) * s * (1-s) },
                             delta = function(t){ rep(0.5,length(t)) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}
tmp <- round(colMeans(results_oneStep,na.rm = T),3)
round((tmp / min(tmp) - 1) * 100,3)
###################################################
###         Operator - ARCH(1)
###################################################
###########
##    1
###########
library(devtools);load_all()
nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'YW'=0,
             'MP90'=0, 'MP95'=0, 'MP99'=0,
             'PW'=0,
             'MLE'=0,
             'Mean'=0, 'Prev'=0, 'Zero'=0)

cat('Fit (',nSims,'): ')
set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
                             alpha = function(s,t,u,v){
                               val_seq <- seq(0,1,length.out=20)
                               K <- 1 #30 #20 #15 #10 #5

                               # tmp <- fda::create.fourier.basis(nbasis = K,rangeval = c(0,1))
                               # basis_vals <- fda::eval.basis(seq(0,1,length.out=length(val_seq)), tmp)
                               #
                               # alp <- matrix(0,nrow=length(val_seq)^2,ncol = length(val_seq)^2)
                               # for(k in 1:K){
                               #   hkk <- basis_vals[,k] %*% t(basis_vals[,k]) /
                               #     nrow(basis_vals)
                               #   # tmp <- max(hkk)
                               #   alp <- alp + kronecker(hkk, hkk) / nrow(hkk)
                               # }
                               # alp <- 225*(alp / K)*20
                               # saveRDS(alp,
                               #         paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'_again3.rds'))
                               alp <- readRDS(
                                 paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'_again3.rds'))
                               alp[(which(s==val_seq)-1)*length(val_seq)+which(t==val_seq),
                                   (which(u==val_seq)-1)*length(val_seq)+which(v==val_seq)]
                             },
                             delta = function(s,t){pmin(t,s) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}
tmp <- round(colMeans(results_oneStep,na.rm = T),3)
round((tmp / min(tmp) - 1) * 100,3)

# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})
# apply(results_oneStep,MARGIN = 2,
#       function(X){
#         x = X[order(-X)]
#         round(mean(x[1:(0.95*length(x))]),3)
#       })
saveRDS(results_oneStep,paste0(getwd(),'/paperData/','opK1.rds'))


###########
##    5
###########
library(devtools);load_all()
nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'YW'=0,
             'MP90'=0, 'MP95'=0, 'MP99'=0,
             'PW'=0,
             'MLE'=0,
             'Mean'=0, 'Prev'=0, 'Zero'=0)

cat('Fit (',nSims,'): ')
set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
                             alpha = function(s,t,u,v){
                               val_seq <- seq(0,1,length.out=20)
                               K <- 5 #30 #20 #15 #10 #5

                               # tmp <- fda::create.fourier.basis(nbasis = K,rangeval = c(0,1))
                               # basis_vals <- fda::eval.basis(seq(0,1,length.out=length(val_seq)), tmp)
                               #
                               # alp <- matrix(0,nrow=length(val_seq)^2,ncol = length(val_seq)^2)
                               # for(k in 1:K){
                               #   hkk <- basis_vals[,k] %*% t(basis_vals[,k]) /
                               #     nrow(basis_vals)
                               #   # tmp <- max(hkk)
                               #   alp <- alp + kronecker(hkk, hkk) / nrow(hkk)
                               # }
                               # alp <- 1200*(alp / K)*20
                               # saveRDS(alp,
                               #         paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'_again.rds'))
                               alp <- readRDS(
                                 paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'_again.rds'))
                               alp[(which(s==val_seq)-1)*length(val_seq)+which(t==val_seq),
                                   (which(u==val_seq)-1)*length(val_seq)+which(v==val_seq)]
                             },
                             delta = function(s,t){pmin(t,s) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}
tmp <- round(colMeans(results_oneStep,na.rm = T),3)
round((tmp / min(tmp) - 1) * 100,3)

# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})
# apply(results_oneStep,MARGIN = 2,
#       function(X){
#         x = X[order(-X)]
#         round(mean(x[1:(0.95*length(x))]),3)
#       })
saveRDS(results_oneStep,paste0(getwd(),'/paperData/','opK5.rds'))

# ###########
# ##    5.2
# ###########
# library(devtools);load_all()
# nSims <- 100
#
# results_oneStep <-
#   data.frame('OR'=rep(NA,nSims),
#              'YW'=NA,
#              'MP'=NA,
#              'PW'=NA,
#              'MLE'=NA,
#              'Mean'=NA,'Prev'=NA,'Zero'=NA)
#
# cat('Fit (',nSims,'): ')
# set.seed(1234567)
# for(i in 1:nSims){
#   cat(paste0(i,', '))
#
#   # non-neg basis function
#   tmp <- suppressMessages(
#     generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
#                              alpha = function(s,t,u,v){
#                                val_seq <- seq(0,1,length.out=20)
#                                K <- 5 #30 #20 #15 #10 #5
#
#                                # tmp <- fda::create.fourier.basis(nbasis = K,rangeval = c(0,1))
#                                # basis_vals <- fda::eval.basis(seq(0,1,length.out=length(val_seq)), tmp)
#                                #
#                                # alp <- matrix(0,nrow=length(val_seq)^2,ncol = length(val_seq)^2)
#                                # for(k in 1:K){
#                                #   hkk <- basis_vals[,k] %*% t(basis_vals[,k]) /
#                                #     nrow(basis_vals)
#                                #   # tmp <- max(hkk)
#                                #   alp <- alp + kronecker(hkk, hkk) / nrow(hkk)
#                                # }
#                                # alp <- 1500*(alp / K)*20
#                                # saveRDS(alp,
#                                #         paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'_again2.rds'))
#                                alp <- readRDS(
#                                  paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'_again2.rds'))
#                                alp[(which(s==val_seq)-1)*length(val_seq)+which(t==val_seq),
#                                    (which(u==val_seq)-1)*length(val_seq)+which(v==val_seq)]
#                              },
#                              delta = function(s,t){pmin(t,s) },
#                              resolution=20)
#   )
#
#
#   results_oneStep[i,] <- tmp
# }
# tmp <- round(colMeans(results_oneStep,na.rm = T),3)
# round((tmp / min(tmp) - 1) * 100,3)
#
# # apply(results_oneStep,MARGIN = 2,
# #       function(X){round(median(X),3)})
# # apply(results_oneStep,MARGIN = 2,
# #       function(X){
# #         x = X[order(-X)]
# #         round(mean(x[1:(0.95*length(x))]),3)
# #       })

###########
##    10
###########
library(devtools);load_all()
nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'YW'=0,
             'MP90'=0, 'MP95'=0, 'MP99'=0,
             'PW'=0,
             'MLE'=0,
             'Mean'=0, 'Prev'=0, 'Zero'=0)


cat('Fit (',nSims,'): ')
set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
                             alpha = function(s,t,u,v){
                               val_seq <- seq(0,1,length.out=20)
                               K <- 10 #30 #20 #15 #10 #5

                               # tmp <- fda::create.fourier.basis(nbasis = K,rangeval = c(0,1))
                               # basis_vals <- fda::eval.basis(seq(0,1,length.out=length(val_seq)), tmp)
                               #
                               # alp <- matrix(0,nrow=length(val_seq)^2,ncol = length(val_seq)^2)
                               # for(k in 1:K){
                               #   hkk <- basis_vals[,k] %*% t(basis_vals[,k]) /
                               #     nrow(basis_vals)
                               #   # tmp <- max(hkk)
                               #   alp <- alp + kronecker(hkk, hkk) / nrow(hkk)
                               # }
                               # alp <- 2500*(alp / K)*20
                               # saveRDS(alp,
                               #         paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'_opt1.rds'))
                               alp <- readRDS(
                                 paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'_opt1.rds'))
                               alp[(which(s==val_seq)-1)*length(val_seq)+which(t==val_seq),
                                   (which(u==val_seq)-1)*length(val_seq)+which(v==val_seq)]
                             },
                             delta = function(s,t){pmin(t,s) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}

tmp <- round(colMeans(results_oneStep,na.rm = T),3)
round((tmp / min(tmp) - 1) * 100,3)

# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})
# apply(results_oneStep,MARGIN = 2,
#       function(X){
#         x = X[order(-X)]
#         round(mean(x[1:(0.95*length(x))]),3)
#       })
saveRDS(results_oneStep,paste0(getwd(),'/paperData/','opK10.rds'))

# ###########
# ##    10.2
# ###########
# library(devtools);load_all()
# nSims <- 100
#
# results_oneStep <-
#   data.frame('OR'=rep(NA,nSims),
#              'OP'=NA,
#              'PW'=NA,
#              'MLE'=NA,
#              'Mean'=NA,'Prev'=NA,'Zero'=NA)
#
#
# cat('Fit (',nSims,'): ')
# set.seed(1234567)
# for(i in 1:nSims){
#   cat(paste0(i,', '))
#
#   # non-neg basis function
#   tmp <- suppressMessages(
#     generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
#                              alpha = function(s,t,u,v){
#                                val_seq <- seq(0,1,length.out=20)
#                                K <- 10 #30 #20 #15 #10 #5
#
#                                # tmp <- fda::create.fourier.basis(nbasis = K,rangeval = c(0,1))
#                                # basis_vals <- fda::eval.basis(seq(0,1,length.out=length(val_seq)), tmp)
#                                #
#                                # alp <- matrix(0,nrow=length(val_seq)^2,ncol = length(val_seq)^2)
#                                # for(k in 1:K){
#                                #   hkk <- basis_vals[,k] %*% t(basis_vals[,k]) /
#                                #     nrow(basis_vals)
#                                #   # tmp <- max(hkk)
#                                #   alp <- alp + kronecker(hkk, hkk) / nrow(hkk)
#                                # }
#                                # alp <- 2750*(alp / K)*20
#                                # saveRDS(alp,
#                                #         paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'_again2.rds'))
#                                alp <- readRDS(
#                                  paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'_again2.rds'))
#                                alp[(which(s==val_seq)-1)*length(val_seq)+which(t==val_seq),
#                                    (which(u==val_seq)-1)*length(val_seq)+which(v==val_seq)]
#                              },
#                              delta = function(s,t){pmin(t,s) },
#                              resolution=20)
#   )
#
#
#   results_oneStep[i,] <- tmp
# }
#
# tmp <- round(colMeans(results_oneStep,na.rm = T),3)
# round((tmp / min(tmp) - 1) * 100,3)
#
# # apply(results_oneStep,MARGIN = 2,
# #       function(X){round(median(X),3)})
# # apply(results_oneStep,MARGIN = 2,
# #       function(X){
# #         x = X[order(-X)]
# #         round(mean(x[1:(0.95*length(x))]),3)
# #       })
# saveRDS(results_oneStep,paste0(getwd(),'/paperData/','opK10_2.rds'))

###########
##    15
###########
###########
##    15
###########
library(devtools);load_all()
nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'YW'=0,
             'MP90'=0, 'MP95'=0, 'MP99'=0,
             'PW'=0,
             'MLE'=0,
             'Mean'=0, 'Prev'=0, 'Zero'=0)

cat('Fit (',nSims,'): ')
set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))
  .predict_ARCH_op_or
  # non-neg basis function
  tmp <- suppressMessages(
    generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
                             alpha = function(s,t,u,v){
                               val_seq <- seq(0,1,length.out=20)
                               K <- 15 #30 #20 #15 #10 #5

                               # tmp <- fda::create.fourier.basis(nbasis = K,rangeval = c(0,1))
                               # basis_vals <- fda::eval.basis(seq(0,1,length.out=length(val_seq)), tmp)
                               #
                               # alp <- matrix(0,nrow=length(val_seq)^2,ncol = length(val_seq)^2)
                               # for(k in 1:K){
                               #   hkk <- basis_vals[,k] %*% t(basis_vals[,k]) /
                               #     nrow(basis_vals)
                               #   # tmp <- max(hkk)
                               #   alp <- alp + kronecker(hkk, hkk) / nrow(hkk)
                               # }
                               # alp <- 3750*(alp / K)*20
                               # saveRDS(alp,
                               #         paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'_opt3.rds'))
                               alp <- readRDS(
                                 paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'_opt3.rds'))
                               alp[(which(s==val_seq)-1)*length(val_seq)+which(t==val_seq),
                                   (which(u==val_seq)-1)*length(val_seq)+which(v==val_seq)]
                             },
                             delta = function(s,t){pmin(t,s) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}
tmp <- round(colMeans(results_oneStep,na.rm = T),3)
round((tmp / min(tmp) - 1) * 100,3)
saveRDS(results_oneStep,paste0(getwd(),'/paperData/','opK15.rds'))

# tmp1 <- (results_oneStep / min(tmp) - 1) * 100
# apply(tmp1,MARGIN = 2,
#       function(X){round(median(X),3)})
# tmp2<- apply(results_oneStep,MARGIN = 2,
#       function(x){
#         x_tmp <- x[order(x)]
#         mean(x_tmp[1:floor(0.99*length(x_tmp))],3)})
# round((tmp2 / min(tmp2) - 1) * 100,3)

###################################################
###         Operator - ARCH(5)
###################################################
###########
##    1
###########
###########
##    5
###########
###########
##    10
###########
###########
##    15
###########
###################################################
###         Operator - GARCH(1,1)
###################################################
###########
##    1
###########
###########
##    5
###########
###########
##    10
###########
###########
##    15
###########
