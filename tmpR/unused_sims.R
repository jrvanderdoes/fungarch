###################################################
###         REPORTED
#       alpha = sum \phi(s)\ten\phi(r) \ten \phi(u)\phi(v)
#       delta = function(s,t){s*(1-s)*t*(1-t) }
###################################################
##########
##    5
##########
library(devtools)
load_all()

nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'OP'=NA,
             'MP'=NA,
             'PW'=NA,
             'MLE'=NA,
             'Mean'=NA,'Prev'=NA,'Zero'=NA)


cat('Fit (',nSims,'): ')
set.seed(12345678) #set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
                             alpha = function(s,t,u,v){
                               const <- 1
                               d <- 5
                               val_seq <- seq(0,1,length.out=20)
                               tmp <- fda::create.fourier.basis(nbasis = d,rangeval = c(0,1))
                               basis_vals <- fda::eval.basis(seq(0,1,length.out=20), tmp)
                               basis_vals_ten <- basis_vals %*% t(basis_vals)/nrow(basis_vals)
                               #alpha <- const * basis_vals_ten %*% t(basis_vals_ten)/nrow(basis_vals_ten)

                               400*basis_vals_ten[which(val_seq==s),which(val_seq==t)] *
                                 basis_vals_ten[which(val_seq==u),which(val_seq==v)]
                             },
                             delta = function(s,t){s*(1-s)*t*(1-t) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}

results_oneStep_old <- results_oneStep
results_oneStep <- na.omit(results_oneStep)

tmp <-round(colMeans(results_oneStep,na.rm = T),3)
# tmp <-
#   round(colMeans(results_oneStep[-c(94,58),],na.rm = T),3)
tmp / min(tmp)
# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})

##########
##    10
##########
library(devtools)
load_all()

nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'OP'=NA,
             'PW'=NA,
             'MLE'=NA,
             'Mean'=NA,'Prev'=NA,'Zero'=NA)


cat('Fit (',nSims,'): ')
set.seed(12345678) #set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
                             alpha = function(s,t,u,v){
                               const <- 1
                               d <- 10
                               val_seq <- seq(0,1,length.out=20)
                               tmp <- fda::create.fourier.basis(nbasis = d,rangeval = c(0,1))
                               basis_vals <- fda::eval.basis(seq(0,1,length.out=20), tmp)
                               basis_vals_ten <- basis_vals %*% t(basis_vals)/nrow(basis_vals)
                               #alpha <- const * basis_vals_ten %*% t(basis_vals_ten)/nrow(basis_vals_ten)

                               400*basis_vals_ten[which(val_seq==s),which(val_seq==t)] *
                                 basis_vals_ten[which(val_seq==u),which(val_seq==v)]
                             },
                             delta = function(s,t){s*(1-s)*t*(1-t) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}

results_oneStep_old <- results_oneStep
results_oneStep <- na.omit(results_oneStep)

tmp <-round(colMeans(results_oneStep,na.rm = T),3)
# tmp <-
#   round(colMeans(results_oneStep[-c(94,58),],na.rm = T),3)
tmp / min(tmp)
# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})

##########
##    15
##########
library(devtools)
load_all()

nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'OP'=NA,
             'PW'=NA,
             'MLE'=NA,
             'Mean'=NA,'Prev'=NA,'Zero'=NA)


cat('Fit (',nSims,'): ')
set.seed(12345678) #set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
                             alpha = function(s,t,u,v){
                               const <- 1
                               d <- 15
                               val_seq <- seq(0,1,length.out=20)
                               tmp <- fda::create.fourier.basis(nbasis = d,rangeval = c(0,1))
                               basis_vals <- fda::eval.basis(seq(0,1,length.out=20), tmp)
                               basis_vals_ten <- basis_vals %*% t(basis_vals)/nrow(basis_vals)
                               #alpha <- const * basis_vals_ten %*% t(basis_vals_ten)/nrow(basis_vals_ten)

                               400*basis_vals_ten[which(val_seq==s),which(val_seq==t)] *
                                 basis_vals_ten[which(val_seq==u),which(val_seq==v)]
                             },
                             delta = function(s,t){s*(1-s)*t*(1-t) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}

results_oneStep_old <- results_oneStep
results_oneStep <- na.omit(results_oneStep)

tmp <-round(colMeans(results_oneStep,na.rm = T),3)
# tmp <-
#   round(colMeans(results_oneStep[-c(94,58),],na.rm = T),3)
tmp / min(tmp)
# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})

##########
##    20/30/40
##########
library(devtools)
load_all()

nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'OP'=NA,
             'PW'=NA,
             'MLE'=NA,
             'Mean'=NA,'Prev'=NA,'Zero'=NA)


cat('Fit (',nSims,'): ')
set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
                             alpha = function(s,t,u,v){
                               val_seq <- seq(0,1,length.out=20)
                               K <- 20

                               # tmp <- fda::create.fourier.basis(nbasis = K,rangeval = c(0,1))
                               # basis_vals <- fda::eval.basis(seq(0,1,length.out=20), tmp)
                               #
                               # alp <- matrix(0,nrow=20^2,ncol = 20^2)
                               # for(k in 1:K){
                               #   hkk <- basis_vals[,k] %*% t(basis_vals[,k]) /
                               #     nrow(basis_vals)
                               #   # tmp <- max(hkk)
                               #   alp <- alp + k*50*kronecker(hkk, hkk) / nrow(hkk)
                               # }
                               # alp <- 100 * alp
                               # saveRDS(alp,
                               #  'C:\\Users\\jerem\\Downloads\\sendGreg\\alp_recheck.rds')
                               alp <- readRDS(
                                       'C:\\Users\\jerem\\Downloads\\sendGreg\\alp_recheck.rds')
                               alp[(which(s==val_seq)-1)*20+which(t==val_seq),
                                   (which(u==val_seq)-1)*20+which(v==val_seq)]
                             },
                             delta = function(s,t){s*(1-s)*t*(1-t) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}

results_oneStep_old <- results_oneStep
results_oneStep <- na.omit(results_oneStep)

# 2.247 2.319 2.382 2.437 2.951 3.029 3.380 (mean)
# 1.252 1.245 1.248 1.242 1.360 1.755 1.769 (median)
round(colMeans(results_oneStep,na.rm = T),3)
tmp<-
  round(colMeans(results_oneStep,na.rm = T),3)
tmp/min(tmp)
# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})

###################################################
###         REPORTED
#       alpha = sum \phi(s)\ten\phi(r) \ten \phi(u)\phi(v)
#       delta = function(s,t){pmin(t,s) }
###################################################
###########
##    5
###########
library(devtools);load_all()
nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'OP'=NA,
             'PW'=NA,
             'MLE'=NA,
             'Mean'=NA,'Prev'=NA,'Zero'=NA)


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
                               # alp <- 1000*(alp / K)*20
                               # saveRDS(alp,
                               #         paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'.rds'))
                               alp <- readRDS(
                                 paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'.rds'))
                               alp[(which(s==val_seq)-1)*length(val_seq)+which(t==val_seq),
                                   (which(u==val_seq)-1)*length(val_seq)+which(v==val_seq)]
                             },
                             delta = function(s,t){pmin(t,s) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}

round(colMeans(results_oneStep,na.rm = T),3)
tmp / min(tmp)
# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})

###########
##    5.
###########
library(devtools);load_all()
nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'OP'=NA,
             'PW'=NA,
             'MLE'=NA,
             'Mean'=NA,'Prev'=NA,'Zero'=NA)


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
                               # alp <- 1250*(alp / K)*20
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

round(colMeans(results_oneStep,na.rm = T),3)
tmp / min(tmp)
# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})


###########
##    10
###########
library(devtools);load_all()
nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'OP'=NA,
             'PW'=NA,
             'MLE'=NA,
             'Mean'=NA,'Prev'=NA,'Zero'=NA)


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
                               # alp <- 2750*(alp / K)*20
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

round(colMeans(results_oneStep,na.rm = T),3)
tmp / min(tmp)
# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})

###########
##    15
###########
library(devtools);load_all()
nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'OP'=NA,
             'PW'=NA,
             'MLE'=NA,
             'Mean'=NA,'Prev'=NA,'Zero'=NA)


cat('Fit (',nSims,'): ')
set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

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
                               # alp <- 3500*(alp / K)*20
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

round(colMeans(results_oneStep,na.rm = T),3)
tmp / min(tmp)
# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})

###########
##    20
###########
library(devtools);load_all()
nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'OP'=NA,
             'PW'=NA,
             'MLE'=NA,
             'Mean'=NA,'Prev'=NA,'Zero'=NA)


cat('Fit (',nSims,'): ')
set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
                             alpha = function(s,t,u,v){
                               val_seq <- seq(0,1,length.out=20)
                               K <- 20 #30 #20 #15 #10 #5

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
                               # alp <- 3500*(alp / K)*20
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

round(colMeans(results_oneStep,na.rm = T),3)
tmp / min(tmp)
# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})

###########
##    30
###########
library(devtools);load_all()
nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'OP'=NA,
             'PW'=NA,
             'MLE'=NA,
             'Mean'=NA,'Prev'=NA,'Zero'=NA)


cat('Fit (',nSims,'): ')
set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
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
                               #   alp <- alp + kronecker(hkk, hkk) / nrow(hkk)
                               # }
                               # alp <- 3750*(alp / K)*20
                               # saveRDS(alp,
                               #         paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'.rds'))
                               alp <- readRDS(
                                 paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'.rds'))
                               alp[(which(s==val_seq)-1)*length(val_seq)+which(t==val_seq),
                                   (which(u==val_seq)-1)*length(val_seq)+which(v==val_seq)]
                             },
                             delta = function(s,t){pmin(t,s) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}

round(colMeans(results_oneStep,na.rm = T),3)
tmp / min(tmp)
# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})


###########
##    40
###########

library(devtools);load_all()
nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'OP'=NA,
             'PW'=NA,
             'MLE'=NA,
             'Mean'=NA,'Prev'=NA,'Zero'=NA)


cat('Fit (',nSims,'): ')
set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitCompareOPARCH(N_fit = 500, N_pred = 250,
                             alpha = function(s,t,u,v){
                               val_seq <- seq(0,1,length.out=20)
                               K <- 40 #30 #20 #15 #10 #5

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
                               # alp <- 2000*(alp / K)*20
                               # saveRDS(alp,
                               #         paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'.rds'))
                               alp <- readRDS(
                                 paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\alp',K,'.rds'))
                               alp[(which(s==val_seq)-1)*length(val_seq)+which(t==val_seq),
                                   (which(u==val_seq)-1)*length(val_seq)+which(v==val_seq)]
                             },
                             delta = function(s,t){pmin(t,s) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}

round(colMeans(results_oneStep,na.rm = T),3)
tmp / min(tmp)
# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})

###################################################
###         REPORTED
#       alpha = sum \phi(s)\ten\phi(r) \ten \phi(u)\phi(v)
#       delta = function(s,t){pmin(t,s) }
###################################################
###########
##    20
###########
###########
library(devtools);load_all()
nSims <- 100

results_oneStep <-
  data.frame('OR'=rep(NA,nSims),
             'OP'=NA,
             'PW'=NA,
             'MLE'=NA,
             'Mean'=NA,'Prev'=NA,'Zero'=NA)


cat('Fit (',nSims,'): ')
set.seed(1234567)
for(i in 1:nSims){
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitCompareOPGARCH(N_fit = 500, N_pred = 250,
                              alpha = function(s,t,u,v){
                                val_seq <- seq(0,1,length.out=20)
                                K <- 20 #30 #20 #15 #10 #5

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
                                # alp <- 3000*(alp / K)*20
                                # saveRDS(alp,
                                #         paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\GARCHalp',K,'_again.rds'))
                                alp <- readRDS(
                                  paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\GARCHalp',K,'_again.rds'))
                                alp[(which(s==val_seq)-1)*length(val_seq)+which(t==val_seq),
                                    (which(u==val_seq)-1)*length(val_seq)+which(v==val_seq)]
                              },
                              beta = function(s,t,u,v){
                                val_seq <- seq(0,1,length.out=20)
                                K <- 20 #30 #20 #15 #10 #5

                                # tmp <- fda::create.fourier.basis(nbasis = K,rangeval = c(0,1))
                                # basis_vals <- fda::eval.basis(seq(0,1,length.out=length(val_seq)), tmp)
                                #
                                # bet <- matrix(0,nrow=length(val_seq)^2,ncol = length(val_seq)^2)
                                # for(k in 1:K){
                                #   hkk <- basis_vals[,k] %*% t(basis_vals[,k]) /
                                #     nrow(basis_vals)
                                #   # tmp <- max(hkk)
                                #   bet <- bet + kronecker(hkk, hkk) / nrow(hkk)
                                # }
                                # bet <- 3000*(bet / K)*20
                                # saveRDS(bet,
                                #         paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\GARCHbet',K,'_again.rds'))
                                bet <- readRDS(
                                  paste0('C:\\Users\\jerem\\Downloads\\sendGreg\\GARCHbet',K,'_again.rds'))
                                alp[(which(s==val_seq)-1)*length(val_seq)+which(t==val_seq),
                                    (which(u==val_seq)-1)*length(val_seq)+which(v==val_seq)]
                              },
                              delta = function(s,t){pmin(t,s) },
                              resolution=20)
  )


  results_oneStep[i,] <- tmp
}
tmp <- round(colMeans(results_oneStep,na.rm = T),3)
tmp / min(tmp)
# apply(results_oneStep,MARGIN = 2,
#       function(X){round(median(X),3)})
