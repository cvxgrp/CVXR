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
                                   args = list(
                                       LinOp(type = 'variable',
                                             size = c(2, 1),
                                             args = list(),
                                             data = 0)
                                   ),
                                   data = NULL)
                         ),
                         data = NULL),
                     constr_id = 1L,
                     size = c(2, 1))
    )


LinOp (id = 2a436cd5-4f9b-4f16-9277-cc5758b8d4c1,
       type = LEQ,
       size = [ 21 ],
       args = LinOp (id = 652ae373-d364-4839-8102-d979c4103c47,
                     type = SUM,
                     size = [ 21 ],
                     args = LinOp (id = b85ae69c-4361-4a54-a5fe-f1ad4408158e,
                                   type = PROMOTE,
                                   size = [ 21 ],
                                   args = LinOp (id = 1fab9af0-56c9-4ffe-ac88-b1246ee10a42,
                                                 type = SCALAR_CONST,
                                                 size = [ 11 ],
                                                 args = ,
                                                 data = [  1  ] )
                                 , data = [    ] )
                     LinOp (id = 58b1ed84-87da-422c-8f6c-e1d7fd21a2c1,
                            type = NEG,
                            size = [ 21 ],
                            args = LinOp (id = a6c3cd62-2064-4ccd-9f9a-8fead749658d,
                                          type = VARIABLE,
                                          size = [ 21 ],
                                          args = ,
                                          data = [  4  ] )
                          , data = [    ] )
                   , data = [    ] )
     , data = [  2  ] )




LinOp (id = f541e27b-3e44-40f2-ba18-2f11015f47e5,
       type = MUL,
       size = [ 11 ],
       args = LinOp (id = 3f82b383-1d80-4f4d-8ea9-c438c67daa62,
                     type = VARIABLE,
                     size = [ 21 ],
                     args = ,
                     data = [  3  ] ),
       data = [  17 19  ] )
