
result_op <- operatorARCH(X = dat_fit, w = 0.01,p = 1 )
# result_op1 <- operatorARCH1(X = dat_fit, w = 0.01,p = 1 )
# result_op2 <- operatorARCH2(X = dat_fit, w = 0.01,p = 1 )
result_op3 <- operatorARCH3(X = dat_fit, w = 0.01,p = 1 )
result_op4 <- operatorARCH4(X = dat_fit, w = 0.01,p = 1 )

estimate_op <-
  .predict_ARCH_op(X2=dat_pred2,
                   alpha=result_op$hat_alpha,
                   delta=result_op$hat_delta,
                   C_eps=C_eps, p=1)
estimate_op3 <-
  .predict_ARCH_op(X2=dat_pred2,
                   alpha=result_op3$hat_alpha,
                   delta=result_op3$hat_delta,
                   C_eps=C_eps, p=1)
estimate_op4 <-
  .predict_ARCH_op(X2=dat_pred2,
                   alpha=result_op4$hat_alpha,
                   delta=result_op4$hat_delta,
                   C_eps=C_eps, p=1)

i <- 3
plot_fd(dat_pred2[,,i])
plot_fd(estimate_op[,,i])
# plot_fd(estimate_op1[,,i])
# plot_fd(estimate_op2[,,i])
plot_fd(estimate_op3[,,i])
plot_fd(estimate_op4[,,i])
