#' Title
#'
#' @param n
#'
#' @return
#' @export
#'
#' @examples
.generate_basis <- function(n){
  C_eps <- .generate_ceps(n)

  # pca_func(C_eps)$rotation
  eigen(C_eps)$vectors*sqrt(n)
}

.generate_ceps <- function(n){
  C_eps <- matrix(NA, n, n)
  times <- 1:n/n
  for (i in 1:n){
    C_eps[i,] <- exp(-abs(times[i]-times)/2)
  }

  C_eps
}


#' Title
#'
#' @param x
#' @param TVE
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
pca_func <- function(x, TVE = 1, ...){
  if(TVE > 1 || TVE < 0) stop('TVE must be in [0,1].',call. = FALSE)

  pc <- stats::prcomp(x=t(x), ...)
  if(TVE==1){
    min_pc <- length(pc$sdev)
  } else{
    min_pc <- min( which(cumsum(pc$sdev^2)/sum(pc$sdev^2) > TVE) )
  }

  ## TODO:: Check out
  new_rot <- pc$rotation[,1:min_pc,drop=FALSE]*sqrt(nrow(x))

  # scores <- matrix(NA, nrow=ncol(x$data),ncol = min_pc)
  # for(i in 1:ncol(x$data)){
  #   for(l in 1:min_pc){
  #     scores[i,l] <- dot_integrate(x$data[,i] %*% t(new_rot[,l]))
  #   }
  # }
  scores <- (t(x-pc$center) %*% new_rot)[,1:min_pc,drop=FALSE]/nrow(x)

  # res = fda::pca.fd(fda::Data2fd(X$data))
  # pc$sdev^2
  # res$values[1:2]
  # View(pc$x); View(res$scores)
  # View(pc$rotation); View(fda::eval.fd(seq(0,1,length.out=30),res$harmonics))

  list(sdev = pc$sdev[1:min_pc]/sqrt(nrow(x)),
       rotation = new_rot,
       center = pc$center,
       scale = pc$scale,
       x = scores,
       orig_pc = pc)
}
