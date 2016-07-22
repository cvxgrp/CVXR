## objective
library(cvxr)

objective <-
    LinOp(type = 'mul',
          size = c(1, 1),
          args = list(
              LinOp(type = 'variable',
                    size = c(2, 1),
                    args = list(),
                    data = 0)
          ),
          data = LinOp(
              type = 'dense_const',
              size = c(1, 2),
              args = list(),
              data = matrix(
                  c(17.,  19.),
                  nrow = 1))
          )

## constraints

constraints <-
    list(
        LinLeqConstr(expr = LinOp(
                         type = 'sum',
                         size = c(2, 1),
                         args = list(
                             LinOp(type = 'promote',
                                   size = c(2, 1),
                                   args = list(
                                       LinOp(type = 'scalar_const',
                                             size = c(1, 1),
                                             args = list(),
                                             data = 1)
                                   ),
                                   data = NULL),
                             LinOp(type = 'neg',
                                   size = c(2, 1),
                                   args = list(LinOp(type = 'variable',
                                                     size = c(2, 1),
                                                     args = list(),
                                                     data = 0)
                                               ),
                                   data = NULL)
                         ), data = NULL),
                     constr_id = 1L,
                     size = c(2, 1))
    )

