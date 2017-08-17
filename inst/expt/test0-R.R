mm <- list(NA,
     list(
         list(
             expr =
                 list(
                     type = 'sum',
                     size = c(2, 1),
                     args = list(
                         list(
                             type = 'promote',
                             size = c(2, 1),
                             args = list(
                                 list(
                                     type = 'scalar_const',
                                     size = c(1, 1),
                                     args = list(),
                                     data = 1),
                                 class = 'LinOp'),
                             data = NULL,
                             class = 'LinOp'),
                         list(
                             type = 'neg',
                             size = c(2, 1),
                             args = list(
                                 list(
                                     type = 'variable',
                                     size = c(2, 1),
                                     args = list(),
                                     data = 5L,
                                     class = 'LinOp')
                             ),
                             data = NULL,
                             class = 'LinOp')),
                     data = NULL,
                     class = 'LinOp'),
             constr_id = 3L,
             size = c(2, 1),
             class = 'LinLeqConstr')
     )
     )
