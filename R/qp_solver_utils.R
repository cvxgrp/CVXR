#' Make ConeDims from constraint list map
#' @param constr_map a list of constraint mappings, typically from a problem
#' @seealso [group_constraints()]
make_cone_dims <- function(constr_map) {
  new("ConeDims",
      zero = Reduce(sum, lapply(constr_map$Zero, size)),
      nonneg = Reduce(sum, lapply(constr_map$NonNeg, size)),
      exp = Reduce(sum, lapply(constr_map$ExpCone, num_cones)),
      soc = unlist(lapply(constr_map$SOC, cone_sizes)),
      psd = unlist(lapply(constr_map$PSD, nrow)),
      p3d = unlist(lapply(constr_map$PowCone3D, function(x) value(x@alpha)))
      )
}

## TODO: The initialize method in dcp2cone.R for ConeDims needs to be replaced by:
setMethod("initialize", "ConeDims",
          function(.Object, zero = 0L, nonneg = 0L, exp = 0L, soc = 0L, psd = 0L, p3d = 0.0)) {
  .Object@zero <- zero
  .Object@nonneg <- nonneg
  .Object@exp <- exp
  .Object@soc <- soc
  .Object@psd <- psd
  .Object@p3d <- p3d
  .Object
})

