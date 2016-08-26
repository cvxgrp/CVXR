intf_size <- function(constant) {
  if(is.numeric(constant) && is.null(dim(constant)) && length(constant) == 1)
    return(c(1,1))
  else if(is.vector(constant))
    return(c(1,length(constant)))
  else if(is.matrix(constant))
    return(dim(constant))
}

intf_sign <- function(constant, tol = 1e-5) {
  c(min(constant) >= -tol, max(constant) <= tol)
}