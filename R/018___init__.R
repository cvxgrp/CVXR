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

# Atom groups.
SOC_ATOMS <- c(
  "GeoMean",
  "Pnorm",
  "QuadForm",
  "QuadOverLin",
  "Power",
  "Huber"
)

EXP_ATOMS <- c(
  "LogSumExp",
  "LogDet",
  "Dcp2Cone.tr_inv_canon",
  "Entr",
  "Exp",
  "KLDiv",
  "Dcp2Cone.rel_entr_canon",
  "Dcp2Cone.von_neumann_entr_canon",
  "Log",
  "Log1p",
  "Logistic",
  "xexp"
)

PSD_ATOMS <- c(
  "LambdaMax",
  "LambdaSumLargest",
  "LogDet",
  "MatrixFrac",
  "NormNuc",
  "SigmaMax"
)

NONPOS_ATOMS <- c(
  "Norm1",
  "Abs",
  "Huber")

## This function is used with several atoms 
apply_with_keepdims <- function(x, fun, axis = NA_real_, keepdims = FALSE) {
  if(is.na(axis))
    result <- fun(x)
  else {
    if(is.vector(x))
      x <- matrix(x, ncol = 1)
    result <- apply(x, axis, fun)
  }

  if(keepdims) {
    new_dim <- dim(x)
    if(is.null(new_dim))
      return(result)
    collapse <- setdiff(1:length(new_dim), axis)
    new_dim[collapse] <- 1
    dim(result) <- new_dim
  }
  result
}

