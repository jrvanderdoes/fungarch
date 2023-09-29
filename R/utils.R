#' Approximate Integral
#'
#' TODO:: UPDATE TO TRAP
#'
#' @param x This can be a vector or matrix
#'
#' @return Numeric estimating intergral value
#'
#' @noRd
.approx_int <- function(x){
  temp_n=NROW(x)
  return((1/temp_n)*sum(x))
}
