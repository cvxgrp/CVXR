list(NA,
     list(
         structure(
             list(
                 expr = structure(
                     list(type = 'sum',
                          size = c(2, 1),
                          args = list(
                              structure(
                                  list(
                                      type = 'promote',
                                      size = c(2, 1),
                                      args = list(
                                          structure(
                                              list(
                                                  type = 'scalar_const',
                                                  size = c(1, 1),
                                                  args = list(
                                                  ),
                                                  data = structure(
                                                      1, .Dim = c(1L,
                                                                  1L)),
                                                  class = 'LinOp'),
                                              .Names = c('type',
                                                         'size',
                                                         'args',
                                                         'data',
                                                         'class'))),
                                      data = NULL,
                                      class = 'LinOp'),
                                  .Names = c('type',
                                             'size',
                                             'args',
                                             'data',
                                             'class')),
                              structure(
                                  list(
                                      type = 'neg',
                                      size = c(2, 1),
                                      args = list(
                                          structure(
                                              list(
                                                  type = 'variable',
                                                  size = c(2, 1),
                                                  args = list(
                                                  ),
                                                  data = 5L,
                                                  class = 'LinOp'),
                                              .Names = c('type',
                                                         'size',
                                                         'args',
                                                         'data',
                                                         'class'))),
                                      data = NULL,
                                      class = 'LinOp'),
                                  .Names = c('type',
                                             'size',
                                             'args',
                                             'data',
                                             'class'))),
                          data = NULL,
                          class = 'LinOp'),
                     .Names = c('type',
                                'size',
                                'args',
                                'data',
                                'class')),
                 constr_id = 3L,
                 size = c(2, 1),
                 class = 'LinLeqConstr'),
             .Names = c('expr',
                        'constr_id',
                        'size',
                        'class'))))



list(
    structure(
        list(
            expr = structure(
                list(
                    type = 'sum_entries',
                    size = c(1, 1),
                    args = list(
                        structure(
                            list(
                                type = 'mul',
                                size = c(2, 1),
                                args = list(
                                    structure(
                                        list(
                                            type = 'variable',
                                            size = c(2, 1),
                                            args = list(
                                            ),
                                            data = 5L,
                                            class = 'LinOp'),
                                        .Names = c('type',
                                                   'size',
                                                   'args',
                                                   'data',
                                                   'class'))),
                                data = structure(
                                    list(

                                        type = 'dense_const',
                                        size = c(2L, 2L),
                                        args = list(
                                        ),
                                        data = structure(
                                            c(1, 0, 0, -1),
                                            .Dim = c(2L, 2L)),
                                        class = 'LinOp'),
                                    .Names = c('type',
                                               'size',
                                               'args',
                                               'data',
                                               'class')),
                                class = 'LinOp'),
                            .Names = c('type',
                                       'size',
                                       'args',
                                       'data',
                                       'class'))),
                    data = NULL,
                    class = 'LinOp'),
                .Names = c('type',
                           'size',
                           'args',
                           'data',
                           'class')),
            constr_id = 7L,
            size = c(1, 1),
            class = 'LinEqConstr'),
        .Names = c('expr',
                   'constr_id',
                   'size',
                   'class')))

[
        LinEqConstr(
            expr=LinOp(type='sum_entries',
                       size=(1, 1),
                       args=[LinOp(type='mul',
                                   size=(2, 1),
                                   args=[LinOp(type='variable',
                                               size=(2, 1),
                                               args=[],
                                               data=0)],
                                   data=LinOp(type='dense_const',
                                              size=(2, 2),
                                              args=[],
                                              data=matrix([[ 1.,  0.],
                                                           [ 0.,  1.]])))],
                       data=None),
            constr_id=3,
            size=(1, 1))]
