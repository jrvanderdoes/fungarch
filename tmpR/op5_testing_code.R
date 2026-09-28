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
for(i in 1:nSims){
  set.seed(1234567+i*123)
  cat(paste0(i,', '))

  # non-neg basis function
  tmp <- suppressMessages(
    generateFitCompareOPARCH5(N_fit = 500, N_pred = 250,

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
                                # alp <- 500*(alp / K)*20
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
