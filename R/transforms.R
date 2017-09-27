linearize <- function(expr) {
  expr <- as.Constant(expr)
  if(is_affine(expr))
    return(expr)
  else {
    tangent <- value(expr)
    if(any(is.na(tangent)))
      stop("Cannot linearize non-affine expression with missing variable values.")
    grad_map <- grad(expr)
    for(var in variables(expr)) {
      grad_var <- grad_map[[as.character(var@id)]]
      if(any(is.na(grad_var)))
        return(NA_real_)
      else if(is_matrix(var)) {
        flattened <- t(Constant(grad_var)) %*% Vec(var - value(var))
        size <- size(expr)
        tangent <- tangent + Reshape(flattened, size[1], size[2])
      } else
        tangent <- tangent + t(Constant(grad_var)) %*% (var - value(var))
    }
  }
  return(tangent)
}
