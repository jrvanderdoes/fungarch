set.seed(1234567)
#set.seed(123)
alpha <- function(t,s){ 16 * t * (1-t) * s * (1-s) }
beta <- function(t,s){ 12 * t * (1-t) * s * (1-s) }
delta <- function(t){ rep(0.1,length(t)) }
dat <- generateFunctionalCondVol(100,1000,'garch',
                                 alpha_par = alpha,
                                 beta_par = beta,
                                 delta_par = delta)

result <- pointwiseGARCH(X = dat$garch_mat, w = 0.05, L=4)

###########################################################

s <- seq(0,1,length.out=nrow(dat$garch_mat))
alpha_true <- matrix(nrow=nrow(dat$garch_mat),ncol=nrow(dat$garch_mat))
beta_true <- matrix(nrow=nrow(dat$garch_mat),ncol=nrow(dat$garch_mat))
for(s1 in 1:length(s)){
  for(s2 in 1:length(s)){
    alpha_true[s1,s2] <- alpha(s[s1],s[s2])
    beta_true[s1,s2] <- beta(s[s1],s[s2])
  }
}
# Compare Alpha
plot3D::persp3D(z=alpha_true,main=parse(text = paste0('~ z == a(s,t)')),
                theta=300,phi=10,shade=0.25,xlab="s",ylab="t",zlab="z")
plot3D::persp3D(z=result$hat_alpha,main=bquote(~z == hat(a)(s,t)),
                theta=300,phi=10,shade=0.25,xlab="s",ylab="t",zlab="z")
# Compare Beta
plot3D::persp3D(z=beta_true,main=parse(text = paste0('~ z == a(s,t)')),
                theta=300,phi=10,shade=0.25,xlab="s",ylab="t",zlab="z")
plot3D::persp3D(z=result$hat_beta,main=bquote(~z == hat(a)(s,t)),
                theta=300,phi=10,shade=0.25,xlab="s",ylab="t",zlab="z")

# Compare Delta
d <- seq(0,1,length.out=nrow(dat$garch_mat))
delta_true <- delta(d)
plot(x=d,delta_true,
     ylim=c(min(delta_true, result$hat_delta),
            max(delta_true, result$hat_delta)),
     type='l')
lines(x=d,y=Re(result$hat_delta), col='red')


.computeGARCHError(dat$garch_mat^2,alpha = alpha_true,
                   beta = beta_true, delta = delta_true)
.computeGARCHError(dat$garch_mat^2,alpha = result$hat_alpha,
                   beta = result$hat_beta, delta = result$hat_delta)
