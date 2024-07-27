## CVXPY SOURCE: none. This is R-specific

#setClassUnion("SolvingChainORNULL", c("SolvingChain", "NULL"))
## We add SolvingChain later to the class union below using setIs function
## see reduction_solvers.R near SolvingChain class defn
setClassUnion("SolvingChainORNULL", c("NULL"))

#setClassUnion("ParamProbORNULL", c("ParamProg", "NULL")) ## THIS SEEMS WRONG, perhaps meant next line
#setClassUnion("ParamProgORNULL", c("ParamQuadProg", "ParamConeProg", "NULL"))
## We add ParamQuadProg, ParamConeProg later to class union below using setIs
## see dcp2cone.R near ParamConeProg class defn and dcp2quad.R near ParamQuadProg
setClassUnion("ParamProgORNULL", c("NULL"))

#setClassUnion("InverseDataORNULL", c("InverseData", "NULL"))
## We add InverseData to InverseDataORNULL later to class union below using setIs
## see reductions.R near InverseData class defn
setClassUnion("InverseDataORNULL", c("NULL"))
