set.seed(123)
dat2 <- generateFunctionalCondVol_ff(resolution = 25,N = 750)

result_pw2 <- pointwiseARCH(X = dat2$garch_mat[,1:500], w = 0.01)
result_op2 <- operatorARCH(X = dat2$garch_mat[,1:500], w = 0.01)

compare_methods_operator(result_pw = result_pw2, result_op = result_op2,
                         result_opor = list('alpha_par'=NULL, 'delta_par'=NULL),
                         X_est1 = dat2$garch_mat[,500:750])

########################

set.seed(1234)
dat3 <- generateFunctionalCondVol_ff(resolution = 21,N = 600)

result_pw3 <- pointwiseARCH(X = dat3$garch_mat[,1:500], w = 0.01)
result_op3 <- operatorARCH(X = dat3$garch_mat[,1:500], w = 0.01)

compare_methods_operator(result_pw = result_pw3, result_op = result_op3,
                         result_opor = list('alpha_par'=NULL, 'delta_par'=NULL),
                         X_est1 = dat3$garch_mat[,500:600])

########################

set.seed(123456)
dat4 <- generateFunctionalCondVol_ff(resolution = 21,N = 1000)

result_pw4 <- pointwiseARCH(X = dat4$garch_mat[,1:750], w = 0.01)
result_op4 <- operatorARCH(X = dat4$garch_mat[,1:750], w = 0.01)

compare_methods_operator(result_pw = result_pw4,
                         result_op = result_op4,
                         result_opor = list('alpha_par'=NULL, 'delta_par'=NULL),
                         X_est1 = dat4$garch_mat[,750:1000])



## predict ahead, out-of-sample, 100 (by-one)
## GARCH PW
##      Curvevol add in
## large alpha
## ADd oracle to Operator
## Check PW in Operator-level
## 500+250 (or so)
## Examples - Setup
## Moving vs Expanding window + horizon size
