nSims <- 1000
ls <- mle <- rep(NA, nSims)

set.seed(123)
for(i in 1:nSims){
  spec <- fGarch::garchSpec(model = list(alpha = c(0.8), beta = 0,delta=4))
  arch1 <- fGarch::garchSim(spec, n = 500)
  result <- acf(arch1^2)
  tmp <- fGarch::garchFit(~garch(1,0),data = arch1)

  ls[i] <- result$acf[2]
  mle[i] <- tmp@fit$matcoef[3,1]
}

library(ggplot2)
ggplot() +
  geom_histogram(aes(x=ls),bins = 40, alpha=0.75 ) +
  geom_vline(aes(xintercept = 0.8), linewidth = 1, color='red', linetype = 'dashed') +
  ylab(NULL) +
  xlab(NULL) +
  theme_bw() +
  theme(axis.text = element_text(size=18),
        axis.text.y = element_blank())

ggplot() +
  geom_histogram(aes(x=mle),bins = 40, alpha=0.75) +
  geom_vline(aes(xintercept = 0.8), linewidth = 1, color='red', linetype = 'dashed') +
  ylab(NULL) +
  xlab(NULL) +
  theme_bw() +
  theme(axis.text = element_text(size=18),
        axis.text.y = element_blank())
