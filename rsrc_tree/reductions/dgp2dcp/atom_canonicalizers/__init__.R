## © Copyright, the CVXR authors

## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at

##     http://www.apache.org/licenses/LICENSE-2.0

## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.

## CVXPY SOURCE: cvxpy/reductions/dgp2dcp/atom_canonicalizers/__init__.py
Dgp2Dcp.CANON_METHODS <- list(AddExpression = Dgp2Dcp.add_canon,
                              Constant = Dgp2Dcp.constant_canon,
                              DivExpression = Dgp2Dcp.div_canon,
                              Exp = Dgp2Dcp.exp_canon,
                              EyeMinusInv = Dgp2Dcp.eye_minus_inv_canon,
                              GeoMean = Dgp2Dcp.geo_mean_canon,
                              GMatMul = Dgp2Dcp.gmatmul_canon,
                              Log = Dgp2Dcp.log_canon,
                              MulExpression = Dgp2Dcp.mulexpression_canon,
                              Multiply = Dgp2Dcp.mul_canon,
                              Norm1 = Dgp2Dcp.norm1_canon,
                              NormInf = Dgp2Dcp.norm_inf_canon,
                              OneMinusPos = Dgp2Dcp.one_minus_pos_canon,
                              PfEigenvalue = Dgp2Dcp.pf_eigenvalue_canon,
                              Pnorm = Dgp2Dcp.pnorm_canon,
                              Power = Dgp2Dcp.power_canon,
                              ProdEntries = Dgp2Dcp.prod_canon,
                              QuadForm = Dgp2Dcp.quad_form_canon,
                              QuadOverLin = Dgp2Dcp.quad_over_lin_canon,
                              Trace = Dgp2Dcp.trace_canon,
                              SumEntries = Dgp2Dcp.sum_canon,
                              XExp = Dgp2Dcp.xexp_canon,
                              Variable = NULL,
                              Parameter = NULL,

                              MaxEntries = EliminatePwl.CANON_METHODS$MaxEntries,
                              MinEntries = EliminatePwl.CANON_METHODS$MinEntries,
                              MaxElemwise = EliminatePwl.CANON_METHODS$MaxElemwise,
                              MinElemwise = EliminatePwl.CANON_METHODS$MinElemwise)

#'
#' DGP canonical methods class.
#'
#' Canonicalization of DGPs is a stateful procedure, hence the need for a class.
#'
#' @rdname DgpCanonMethods-class
.DgpCanonMethods <- setClass("DgpCanonMethods", representation(.variables = "list", .parameters = "list"), prototype(.variables = list(), .parameters = list()), contains = "list")
DgpCanonMethods <- function(...) { .DgpCanonMethods(...) }

#' @param x A \linkS4class{DgpCanonMethods} object.
#' @describeIn DgpCanonMethods Returns the name of all the canonicalization methods
setMethod("names", signature(x = "DgpCanonMethods"), function(x) { names(Dgp2Dcp.CANON_METHODS) })

setMethod("match", signature(x = "character", table = "DgpCanonMethods"), function(x, table, nomatch = NA_integer_, incomparables = NULL) {
  return(x %in% names(Dgp2Dcp.CANON_METHODS))
})

## # TODO: How to implement this with S4 setMethod? Signature is x = "DgpCanonMethods", i = "character", j = "missing".
## '[[.DgpCanonMethods' <- function(x, i, j, ..., exact = TRUE) { do.call("$", list(x, i)) }

setMethod(f = "[[", signature(x = "DgpCanonMethods", i = "character", j = "missing"),
          definition = function(x, i, j, ..., exact = TRUE) {
              # Ensure index 'i' is of character type
              if (!is.character(i)) {
                  stop("index 'i' must be of type character")
              }
              
              # Attempt to fetch the value from the .variables slot
              if (i %in% names(x@.variables)) {
                  return(x@.variables[[i]])
              }
              
              # If not found, attempt to fetch the value from the .parameters slot
              if (i %in% names(x@.parameters)) {
                  return(x@.parameters[[i]])
              }
              
              # If the index is not found in either slots, return NULL or throw an error
              stop(sprintf("Variable or parameter '%s' not found in 'DgpCanonMethods' object", i))
          })


#' @param name The name of the atom or expression to canonicalize.
#' @describeIn DgpCanonMethods Returns either a canonicalized variable or
#'  a corresponding Dgp2Dcp canonicalization method
setMethod("$", signature(x = "DgpCanonMethods"), function(x, name) {
  if(name == "Variable")
    return(function(...) { variable_canon(x, ...) })
  else if(name == "Parameter")
    return(function(...) { parameter_canon(x, ...) })
  else
    return(Dgp2Dcp.CANON_METHODS[[name]])
})

setMethod("variable_canon", "DgpCanonMethods", function(object, variable, args) {
  # args <- NULL
  # Swaps out positive variables for unconstrained variables.
  vid <- id(variable)
  vid_char <- as.character(vid)
  if(vid_char %in% names(object@.variables))
    return(list(object, object@.variables[[vid_char]], list()))
  else {
    log_variable <- new("Variable", dim = dim(variable), var_id = vid)
    object@.variables[[vid_char]] <- log_variable
    return(list(object, log_variable, list()))
  }
})

setMethod("parameter_canon", "DgpCanonMethods", function(object, parameter, args) {
  # args <- NULL
  # Swaps out positive parameters for unconstrained variables.
  pid <- id(parameter)
  pid_char <- as.character(pid)
  if(pid_char %in% names(object@.parameters))
    return(list(object, object@.parameters[[pid_char]], list()))
  else {
    log_parameter <- new("Parameter", dim = dim(parameter), name = name(parameter), value = base::log(value(parameter)))
    object@.parameters[[pid_char]] <- log_parameter
    return(list(object, log_parameter, list()))
  }
})
