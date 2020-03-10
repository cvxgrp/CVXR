library(testthat)
library(CVXR)

test_check("CVXR", filter="^constant_atoms")
