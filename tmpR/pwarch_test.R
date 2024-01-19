# NOrm on alpha <1, delta>0
# set delta to 0 if below
# Try trainign and test in fitting
# Try m2 as mean (X^2)

set.seed(123)
alpha <- function(t,s){ 16 * t * (1-t) * s * (1-s) }
delta <- function(t){ rep(0.1,length(t)) }
dat <- generateFunctionalCondVol(resolution = 25,N = 600,'arch',
                                 alpha_par = alpha,
                                 delta_par = delta)

result_pw <- pointwiseARCH(X = dat$garch_mat[,1:500], w = 0.01)
#result_op <- operatorARCH(X = dat$garch_mat, w = 0.01)

compare_methods_pointwise(result_pw, get_oracle(dat), dat$garch_mat[,500:600])

#################################################

set.seed(1234567)
alpha1 <- function(t,s){ 16 * t * (1-t) * s * (1-s) }
delta1 <- function(t){ rep(0.01,length(t)) }
dat1 <- generateFunctionalCondVol(resolution = 25,N = 1250,'arch',
                                 alpha_par = alpha1,
                                 delta_par = delta1)

result_pw1 <- pointwiseARCH(X = dat1$garch_mat[,1:1000], w = 0.01)
#result_op1 <- operatorARCH(X = dat1$garch_mat[,1:1000], w = 0.01)

compare_methods_pointwise(result_pw1, get_oracle(dat1),
                          dat1$garch_mat[,1000:1250])

#################################################

# Show Seeds
# Ask about big c on alpha
set.seed(1234) #set.seed(12345)
alpha0 <- function(t,s){ 32 * t * (1-t) * s * (1-s) }
delta0 <- function(t){ rep(0.5,length(t)) }
dat0 <- generateFunctionalCondVol(resolution = 20,N = 750,'arch',
                                  alpha_par = alpha0,
                                  delta_par = delta0)

result_pw0 <- pointwiseARCH(X = dat0$garch_mat[,1:500], w = 0.01) # better than w=0.01
#result_pw0 <- pointwiseARCH(X = dat0$garch_mat[,1:500], w = 0.07)
#result_op1 <- operatorARCH(X = dat1$garch_mat[,1:1000], w = 0.01)

compare_methods_pointwise(result_pw0, get_oracle(dat0), dat0$garch_mat[,500:750])

#################################################

# Show Seeds
# Ask about big c on alpha
set.seed(1234567) #set.seed(12345)
alpha00 <- function(t,s){ 32 * t * (1-t) * s * (1-s) }
delta00 <- function(t){ rep(0.5,length(t)) }
dat00 <- generateFunctionalCondVol(resolution = 20,N = 6000,'arch',
                                  alpha_par = alpha00,
                                  delta_par = delta00)

result_pw00 <- pointwiseARCH(X = dat00$garch_mat[,1:5000], w = 0.05)
#result_op1 <- operatorARCH(X = dat1$garch_mat[,1:1000], w = 0.01)

# compare_methods_pointwise(result_pw00, get_oracle(dat00),
#                           dat00$garch_mat[,4000:5000])
# compare_methods_pointwise(result_pw_out00, get_oracle(dat00),
#                           dat00$garch_mat[,4000:5000])
compare_methods_pointwise(result_pw00, get_oracle(dat00),
                          dat00$garch_mat[,5000:6000])

###########################################################
result <- result_pw00
dat <- dat00

s <- seq(0,1,length.out=nrow(dat$garch_mat))
alpha_true <- matrix(nrow=nrow(dat$garch_mat),ncol=nrow(dat$garch_mat))
for(s1 in 1:length(s)){
  for(s2 in 1:length(s)){
    alpha_true[s1,s2] <- dat$params$alpha_par(s[s1],s[s2])
  }
}
# Compare Alpha
plot3D::persp3D(z=alpha_true,main=parse(text = paste0('~ z == a(s,t)')),
                theta=300,phi=10,shade=0.25,xlab="s",ylab="t",zlab="z")
plot3D::persp3D(z=result$hat_alpha,main=bquote(~z == hat(a)(s,t)),
                theta=300,phi=10,shade=0.25,xlab="s",ylab="t",zlab="z")

# Compare Delta
d <- seq(0,1,length.out=nrow(dat$garch_mat))
delta_true <- dat$params$delta_par(d)
plot(x=d,delta_true,
     ylim=c(min(delta_true, result$hat_delta),
            max(delta_true, result$hat_delta)),
     type='l')
lines(x=d,y=as.numeric(result$hat_delta), col='red')
#lines(x=d,y=as.numeric(rowMeans(result$hat_delta)), col='red')
#lines(x=seq(0,1,length.out=1000),y=as.numeric(result$hat_delta), col='red')


# Compare Delta II
d <- seq(0,1,length.out=nrow(dat$garch_mat))
delta_true <- sapply(1:nrow(dat$garch_mat),function(x,d){dat$params$delta_par(d)},d=d)
plot3D::persp3D(z=delta_true,main=parse(text = paste0('~ z == a(s,t)')),
                theta=300,phi=10,shade=0.25,xlab="s",ylab="t",zlab="z")
plot3D::persp3D(z=result$hat_delta,main=bquote(~z == hat(a)(s,t)),
                theta=300,phi=10,shade=0.25,xlab="s",ylab="t",zlab="z")

