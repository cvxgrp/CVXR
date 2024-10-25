## These are class unions we need for our own classes to extend.  The
## need for these unions arises in specific code files, but by
## defining all the needed ones first here in the simplest form, we
## can ensure that they are available for code using these classes as
## signatures in method definitions.  THEN, as the various classes are
## defined, we can add those classes to the union by using `setIs`
## (see ?setIs in R.)

## These class unions are based on either base R types or on packages
## imported. So they can always go first in any collation order.

## BEGIN GROUP
setClassUnion("NumORNULL", c("numeric", "NULL"))

setClassUnion("NumORLogical", c("logical", "numeric"))

setClassUnion("S4ORNULL", c("S4", "NULL"))

setClassUnion("ListORNULL", c("list", "NULL"))

setClassUnion("NumORgmp", c("numeric", "bigq", "bigz"))

setClassUnion("ConstSparseVal", c("CsparseMatrix", "TsparseMatrix"))

setClassUnion("ConstVal", c("ConstSparseVal", "data.frame", "matrix", "numeric",
                            "complex", "dMatrix", "bigq", "bigz"))

setClassUnion("ConstValORNULL", c("ConstVal", "NULL"))

## END GROUP

## BEGIN GROUP
## This needed in constraints/constraint.R and does not get used anywhere else
## We begin by just adding plain list and add constraint class later after it is defined
setClassUnion("ListORConstr", c("list"))
## END GROUP 

## BEGIN GROUP
## This is needed in expressions/expression.R and does not get used anywhere else
setClassUnion("ConstValORExpr", "ConstVal") ## Begin with just ConstVal
setClassUnion("ConstValListORExpr", c("ConstVal", "list")) ## Begin with what's available here
setClassUnion("ListORExpr", "list") ## Begin with plain list
## END GROUP

## Begin GROUP
## This is needed in reductions/solvers/solving_chain.R and problems/problem.R
setClassUnion("SolvingChainORNULL", "NULL")
## END GROUP

## Begin GROUP
## We add ParamQuadProg, ParamConeProg later to class union below using setIs
## see dcp2cone.R near ParamConeProg class defn and dcp2quad.R near ParamQuadProg
setClassUnion("ParamProgORNULL", "NULL")  # Begin with NULL
## END GROUP

## Begin GROUP
## We add InverseData to InverseDataORNULL later to class union below using setIs
## see reductions.R near InverseData class defn
setClassUnion("InverseDataORNULL", "NULL")
## End GROUP


## Begin GROUP
## This is updated in problems/problem.R
setClassUnion("SolverStatsORNULL", "NULL")
setClassUnion("SizeMetricsORNULL", "NULL")
## End GROUP

## Begin GROUP
## This is updated in problems/problem.R
setClassUnion("SolutionORList", "list")
setClassUnion("ProblemORNULL", "NULL")
## End GROUP

## Begin GROUP
## This is updated in reductions/utilities.R
setClassUnion("ReducedMatORNULL", "NULL")
## End GROUP

## Begin GROUP
## This is further updated in reductions/solvers/solver.R
setClassUnion("ReductionSolverORNULL", "NULL")
## End GROUP







