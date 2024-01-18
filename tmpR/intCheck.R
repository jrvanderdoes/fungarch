int <- 0
s_vals <- seq(0,1,0.00001)
for(s in s_vals){
  int <- int + integrate(function(t,s){ 36* t * (1-t) * s * (1-s) },lower = 0,upper = 1,s=s)[[1]]
}
int / length(s_vals)
1/(int / length(s_vals))
