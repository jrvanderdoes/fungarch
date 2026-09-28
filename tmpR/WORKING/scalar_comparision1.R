library(devtools);load_all()

tmp <- .generate_ceps(20)
dat <- generateFunctionalARCH_kl(
  delta = 10 *(tmp[,1] %o% tmp[,1]),
  resolution = 20,
  N = 1000,
  burnin = 100,
  p = 1, K=1,scale = 0.8)
X = dat$data1
w = 0.001
p = 1
num_comp = 3


plot( as.vector( t(X) %*% tmp[,1] ), type='l')

acf((t(X) %*% tmp[,1])^2)

spec <- fGarch::garchSpec(model = list(alpha = c(0.8), beta = 0,delta=4))
arch1 <- fGarch::garchSim(spec, n = 1000)
result <- acf(arch1^2)
tmp <- fGarch::garchFit(~garch(1,0),data = arch1)

result$acf[2]
tmp@fit$matcoef[3,1]
