intf_size <- function(constant) {
  if((is.null(dim(constant)) && length(constant) == 1) || all(dim(constant) == c(1,1)))
    return(c(1,1))
  else if(is.vector(constant))
    return(c(1,length(constant)))
  else if(is.matrix(constant))
    return(dim(constant))
}

intf_sign <- function(constant, tol = 1e-5) {
  c(min(constant) >= -tol, max(constant) <= tol)
}

intf_scalar_value <- function(constant) {
  if(is.list(constant))
    constant[[1]]
  else
    as.numeric(constant)
}