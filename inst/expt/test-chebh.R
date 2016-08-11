##### THIS EXAMPLE WORKS!!!!

##The main obj
library(cvxr)

C_objective <- CVXcanon.LinOp$new(
    type = 'NEG',
    size = c(1, 1),
    args = list(
        CVXcanon.LinOp$new(type = 'NEG',
                           size = c(1,1),
                           args = list(CVXcanon.LinOp$new(type = 'VARIABLE',
                                                          size = c(1, 1),
                                                          data = matrix(0, 1))
                                       )
                           )
    )
)

C_constraints <- CVXcanon.LinOpVector$new()

c1 <- CVXcanon.LinOp$new(type = 'LEQ',
                         size = c(1, 1),
                         args= list(
                             CVXcanon.LinOp$new(type = 'SUM',
                                                size = c(1, 1),
                                                args = list(
                                                    CVXcanon.LinOp$new(type = 'MUL',
                                                                       size = c(1, 1),
                                                                       args = list(
                                                                           CVXcanon.LinOp$new(type = 'VARIABLE',
                                                                                              size = c(2, 1),
                                                                                              data = matrix(1, 1))
                                                                       ),
                                                                       data = matrix(c(2, 1), nrow = 1)),
                                                    CVXcanon.LinOp$new(type = 'MUL',
                                                                       size = c(1, 1),
                                                                       args = list(
                                                                           CVXcanon.LinOp$new(type = 'VARIABLE',
                                                                                              size = c(1, 1),
                                                                                              data = matrix(0, 1))
                                                                       ),
                                                                       data = matrix(2.23606797749979, 1)),
                                                    CVXcanon.LinOp$new(type = 'SCALAR_CONST',
                                                                       size = c(1, 1),
                                                                       data = matrix(-1, 1))
                                                )
                                                )
                         ),
                         data = matrix(2, 1) ## constr_id = 2
                         )

c2 <- CVXcanon.LinOp$new(type = 'LEQ',
                         size = c(1, 1),
                         args= list(
                             CVXcanon.LinOp$new(type = 'SUM',
                                                size = c(1, 1),
                                                args = list(
                                                    CVXcanon.LinOp$new(type = 'MUL',
                                                                       size = c(1, 1),
                                                                       args = list(
                                                                           CVXcanon.LinOp$new(type = 'VARIABLE',
                                                                                              size = c(2, 1),
                                                                                              data = matrix(1, 1))
                                                                       ),
                                                                       data = matrix(c(2, -1), nrow = 1)),
                                                    CVXcanon.LinOp$new(type = 'MUL',
                                                                       size = c(1, 1),
                                                                       args = list(
                                                                           CVXcanon.LinOp$new(type = 'VARIABLE',
                                                                                              size = c(1, 1),
                                                                                              data = matrix(0, 1))
                                                                       ),
                                                                       data = matrix(2.23606797749979, 1)),
                                                    CVXcanon.LinOp$new(type = 'SCALAR_CONST',
                                                                       size = c(1, 1),
                                                                       data = matrix(-1, 1))
                                                )
                                                )
                         ),
                         data = matrix(3, 1) ## constr_id = 3
                         )

c3 <- CVXcanon.LinOp$new(type = 'LEQ',
                         size = c(1, 1),
                         args= list(
                             CVXcanon.LinOp$new(type = 'SUM',
                                                size = c(1, 1),
                                                args = list(
                                                    CVXcanon.LinOp$new(type = 'MUL',
                                                                       size = c(1, 1),
                                                                       args = list(
                                                                           CVXcanon.LinOp$new(type = 'VARIABLE',
                                                                                              size = c(2, 1),
                                                                                              data = matrix(1, 1))
                                                                       ),
                                                                       data = matrix(c(-1, 2), nrow = 1)),
                                                    CVXcanon.LinOp$new(type = 'MUL',
                                                                       size = c(1, 1),
                                                                       args = list(
                                                                           CVXcanon.LinOp$new(type = 'VARIABLE',
                                                                                              size = c(1, 1),
                                                                                              data = matrix(0, 1))
                                                                       ),
                                                                       data = matrix(2.23606797749979, 1)),
                                                    CVXcanon.LinOp$new(type = 'SCALAR_CONST',
                                                                       size = c(1, 1),
                                                                       data = matrix(-1, 1))
                                                )
                                                )
                         ),
                         data = matrix(4, 1) ## constr_id = 4
                         )


c4 <- CVXcanon.LinOp$new(type = 'LEQ',
                         size = c(1, 1),
                         args= list(
                             CVXcanon.LinOp$new(type = 'SUM',
                                                size = c(1, 1),
                                                args = list(
                                                    CVXcanon.LinOp$new(type = 'MUL',
                                                                       size = c(1, 1),
                                                                       args = list(
                                                                           CVXcanon.LinOp$new(type = 'VARIABLE',
                                                                                              size = c(2, 1),
                                                                                              data = matrix(1, 1))
                                                                       ),
                                                                       data = matrix(c(-1, -2), nrow = 1)),
                                                    CVXcanon.LinOp$new(type = 'MUL',
                                                                       size = c(1, 1),
                                                                       args = list(
                                                                           CVXcanon.LinOp$new(type = 'VARIABLE',
                                                                                              size = c(1, 1),
                                                                                              data = matrix(0, 1))
                                                                       ),
                                                                       data = matrix(2.23606797749979, 1)),
                                                    CVXcanon.LinOp$new(type = 'SCALAR_CONST',
                                                                       size = c(1, 1),
                                                                       data = matrix(-1, 1))
                                                )
                                                )
                         ),
                         data = matrix(5, 1) ## constr_id = 5
                         )

C_constraints$push_back(c1)
C_constraints$push_back(c2)
C_constraints$push_back(c3)
C_constraints$push_back(c4)

cat("Objective:\n")
print(C_objective$toString())

cat("Constraints:\n")
print(C_constraints$toString())

C_sense <- 1L ## CVXcanon.MAXIMIZE

C_opts <- list(verbose = 1L)
cvxCanon <- CVXcanon$new()
solution <- cvxCanon$solve(C_sense, C_objective, C_constraints, C_opts)
print(solution$optimal_value * -1)
print(solution$primal_values)

###


