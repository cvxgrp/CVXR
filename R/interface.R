intf_size <- function(constant) {
  if((is.null(dim(constant)) && length(constant) == 1) || 
     (!is.null(dim(constant)) && all(dim(constant) == c(1,1))))
    return(c(1,1))
  else if(is.atomic(constant))
    return(c(1,length(constant)))
  else if(is.matrix(constant) || is(constant, "Matrix"))
    return(dim(constant))
  else
    stop("Unknown class: ", class(constant))
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

intf_convert_if_scalar <- function(constant) {
  if(all(intf_size(constant) == c(1,1)))
    intf_scalar_value(constant)
  else
    constant
}