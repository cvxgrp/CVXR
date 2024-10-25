## CVXPY SOURCE: cvxpy/reductions/cone2cone/.exotic2common.py

# An "exotic" cone is defined as any cone that isn't supported by ParamConeProg.
# If ParamConeProg is updated to support more cones, then it may be necessary to change this file.
EXOTIC_CONES <- list("PowConeND" = list("PowCone3D"))

# con : PowConeND. We can extract metadata from this, e.g., con@alpha and con@axis.
# args : tuple of length two with W = args[[1]], z = args[[2]].
pow_nd_canon <- function(con, args) {
  data <- get_data(con)
  alpha <- data[[1]]
  axis <- data[[2]]
  alpha <- value(alpha)
  W <- args[[1]]
  z <- args[[2]]

  if(axis == 1) {
    W <- t(W)
    alpha <- t(alpha)
  }

  if(ndim(W) == 1) {
    W <- Reshape(W, c(size(W), 1))
    alpha <- Reshape(alpha, c(size(W), 1))
  }

  W_dim <- dim(W)
  n <- W_dim[1]
  k <- W_dim[2]
  if(n == 2)
    can_canon <- PowCone3D(W[1,], W[2,], z, alpha[1,])
  else {
    Tvar <- Variable(n-2, k)

    x_3d <- list()
    y_3d <- list()
    z_3d <- list()
    alpha_3d <- list()
    for(j in seq_len(k)) {
      x_3d <- c(x_3d, W[seq_len(n-1),j])
      y_3d <- c(y_3d, Tvar[,j])
      y_3d <- c(y_3d, W[n,j])
      z_3d <- c(z_3d, z[j])
      z_3d <- c(z_3d, Tvar[,j])

      r_nums <- alpha[,j]
      r_dens <- rev(base::cumsum(rev(r_nums)))   # Equivalent to sapply(seq(n), function(i) { sum(alpha[seq(i,n), j]) }).
      r <- r_nums / r_dens
      alpha_3d <- c(alpha_3d, r[seq_len(n-1)])
    }

    # TODO: Ideally, we should construct x, y, z, alpha_p3d by applying suitable sparse matrices to W, z, Tvar,
    # rather than using the HStack atom. (HStack will probably result in longer compile times).
    x_3d <- do.call(HStack, x_3d)
    y_3d <- do.call(HStack, y_3d)
    z_3d <- do.call(HStack, z_3d)
    alpha_p3d <- do.call(HStack, alpha_3d)

    can_con <- PowCone3D(x_3d, y_3d, z_3d, alpha_p3d)
  }

  # Return a single PowCone3D constraint defined over all auxiliary variables needed for the reduction to go through. There are no "auxiliary constraints" beyond this one.
  return(list(can_con, list()))
}

#'
#' The Exotic2Common class.
#'
#' This class represents a reduction of an exotic cone to a common cone constraint.
#'
#' @rdname Exotic2Common-class
.Exotic2Common <- setClass("Exotic2Common", contains = "Canonicalization")

Exotic2Common <- function(problem = NULL) { .Exotic2Common(problem = problem) }

Exotic2Common.CANON_METHODS <- list("PowConeND" = pow_nd_canon)

setMethod("initialize", "Exotic2Common", function(.Object, ...) {
  callNextMethod(.Object, ..., canon_methods = Exotic2Common.CANON_METHODS)
})
