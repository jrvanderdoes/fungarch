#' Approximate Integral
#'
#' TODO:: UPDATE TO TRAP
#'
#' @param x This can be a vector or matrix
#'
#' @return
#' @export
#'
#' @examples
.approx_int <- function(x){
  temp_n=NROW(x)
  return((1/temp_n)*sum(x))
}
