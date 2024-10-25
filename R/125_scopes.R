## cvxpy note: cvxpy/utilities/scopes.py

.CVXR_options$dpp_scope_active <- FALSE

## Context manager for DPP curvature analysis
## When this scope is active, parameters are affine, not constant. For example, if `param` is a Parameter, then
##  if `dpp_scope()` is TRUE:
## `cat("param is constant: ", param.is_constant())`
## `cat("param is affine: ", param.is_affine())`
##  would print
##  `param is constant: FALSE`
##  `param is affine: TRUE`
##
## def dpp_scope() -> Generator[None, None, None]:
##     global _dpp_scope_active
##     prev_state = _dpp_scope_active
##     _dpp_scope_active = True
##     yield
##     _dpp_scope_active = prev_state


#' Returns TRUE if a `dpp_scope` is active
#' @return TRUE OR FALSE
dpp_scope_active <- function() {
  .CVXR_options$dpp_scope_active
}
