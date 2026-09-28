operator <- function(X, Y){
  X %*% as.vector(Y) / nrow(X)
}

operator_vech <- function(X, Y){
  X %*% ks::vech(Y) / nrow(X)
}
# tensor <- function(X,Y){
#   X %o% Y
# }

tensor_vec <- function(X,Y){
  t( as.vector(X) %o% as.vector(Y) )
  # kronecker(as.vector(X),t(as.vector(Y)))
}
tensor_vech <- function(X,Y){
  t( ks::vech(X) %o% ks::vech(Y) )
}

adjacent <- function(X,Y){
  # ifelse(!is.null(dim(X)),
  #        X %*% Y / nrow(X),
  #        X %*% Y / length(X) )
  X %*% Y / ncol(X)
}
