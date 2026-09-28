

my_int <- function(times,prev_val){
  dot_integrate(times*(1-times)*prev_val)
}

TMP <- sapply(times,delta_par,t=times) +
  ( times*(1-times) ) %*% t( times*(1-times) ) *
   my_int(times, sim_garch_matrix[,i-1]) *
   my_int(times, sim_garch_matrix[,i-1])
#View(TMP)


prev_tensor <- sim_garch_matrix[,i-1] %*% t(sim_garch_matrix[,i-1])
for(j in 1:resolution){
  # Over all t
  for(k in 1:resolution){
    vector_op[,,k] <- vector_list[[j]][[k]] %*% prev_tensor / length(times)
  }
  # TODO:: Convert to trap int
  vector_op1[,,j] <- rowMeans(vector_op, dims=2)
  # vector_op1_tmp[,,j] <- (vector_op[,,1] + vector_op[,,dim(vector_op)[3]] +
  #           rowSums(2*vector_op[,,2:(dim(vector_op)[3]-1)],dims=2) ) /
  #   (2*dim(vector_op)[3])
}
# Formula (2)
#  delta + sum alpha(X2) + sum beta(Sigma2)
# TODO:: Convert to trap int
sim_sigma2_matrix[,,i] <-
  sapply(times,delta_par,t=times) + rowMeans(vector_op1, dims=2)


View(sim_sigma2_matrix[,,i])
View(TMP)
