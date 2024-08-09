## These are class unions we need for our own classes to extend

## These class unions are based on either base R types or on packages imported. So they can always go first
## in the collation order
## BEGIN GROUP
setClassUnion("NumORNULL", c("numeric", "NULL"))

setClassUnion("NumORLogical", c("logical", "numeric"))

setClassUnion("S4ORNULL", c("S4", "NULL"))

setClassUnion("ListORNULL", c("list", "NULL"))

setClassUnion("NumORgmp", c("numeric", "bigq", "bigz"))

setClassUnion("ConstSparseVal", c("CsparseMatrix", "TsparseMatrix"))

setClassUnion("ConstVal", c("ConstSparseVal", "data.frame", "matrix", "numeric", "complex", "dMatrix", "bigq", "bigz"))
## END GROUP

## BEGIN GROUP
## constraints/constraint.R
setClassUnion("ListORConstr", c("list", "Constraint"))
## END GROUP 

## BEGIN GROUP
## 
setClassUnion("ConstValORExpr", c("ConstVal", "Expression"))
setClassUnion("ConstValORNULL", c("ConstVal", "NULL"))
setClassUnion("ConstValListORExpr", c("ConstVal", "list", "Expression"))
setClassUnion("ListORExpr", c("list", "Expression"))
## END GROUP
#setClassUnion("SolvingChainORNULL", c("SolvingChain", "NULL"))

setClassUnion("SolvingChainORNULL", c("NULL"))

#setClassUnion("ParamProbORNULL", c("ParamProg", "NULL")) ## THIS SEEMS WRONG, perhaps meant next line

#setClassUnion("ParamProgORNULL", c("ParamQuadProg", "ParamConeProg", "NULL"))

setClassUnion("ParamProgORNULL", c("NULL"))

#setClassUnion("InverseDataORNULL", c("InverseData", "NULL"))

setClassUnion("InverseDataORNULL", c("NULL"))

setClassUnion("SolverStatsORNULL", c("SolverStats", "NULL"))

setClassUnion("SizeMetricsORNULL", c("SizeMetrics", "NULL"))

setClassUnion("SolutionORList", c("Solution", "list"))

setClassUnion("ReductionSolverORNULL", c("ReductionSolver", "NULL"))

setClassUnion("ReducedMatORNULL", c("ReducedMat", "NULL"))

setClassUnion("ProblemORNULL", c("Problem", "NULL"))
