## CVXPY SOURCE: cvxpy/reductions/dgp2dcp/util.py

# TODO: Do we need this in R? The Python sum function is a reduction with initial value 0.0, resulting in a non-DGP expression.
Dgp2Dcp.explicit_sum <- function(expr) {
  x <- Vec(expr)
  summation <- x[1]
  x_len <- size(x)
  if(x_len > 1) {
    for(i in seq(2, x_len))
      summation <- summation + x[i]
  }
  return(summation)
}

