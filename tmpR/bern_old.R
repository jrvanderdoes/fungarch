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
                             alpha = function(s,t,u,v,const = 7000){
                               d <- 25
                               bern <- function (k, n, x)
                               {
                                 stopifnot(is.numeric(k), is.numeric(n), is.numeric(x))
                                 if (floor(n) != ceiling(n) || n < 0)
                                   stop("Argument 'n' must be a nonnegative integer.")
                                 if (floor(k) != ceiling(k) || k < 0 || k > n)
                                   stop("Argument 'k' must be an integer between 0 and 'n'.")
                                 if (any(x < -1 | x > 1))
                                   stop("Elements of argument 'x' must all lie in [0, 1].")
                                 choose(n, k) * x^k * (1 - x)^(n - k)
                               }

                               tmp <- 0
                               for(i in 1:d){
                                 tmp <- tmp +
                                   bern(i,d,s) * bern(i,d,t) *
                                   bern(i,d,u) * bern(i,d,v)
                               }
                               const * tmp
                             },
                             delta = function(s,t){s*(1-s)*t*(1-t) },
                             resolution=20)
  )


  results_oneStep[i,] <- tmp
}
# saveRDS(results_oneStep, 'C://Users//jerem//Downloads//op_onestep.rds')
results_oneStep_old <- results_oneStep
results_oneStep <- na.omit(results_oneStep)

round(colMeans(results_oneStep,na.rm = T),3)
apply(results_oneStep,MARGIN = 2,
      function(X){round(median(X),3)})
