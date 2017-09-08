[LinEqConstr(expr=LinOp(type='sum', size=(2, 2),
                        args=[LinOp(type='index', size=(2, 2),
                                    args=[LinOp(type='variable', size=(3, 3), args=[], data=1)],
                                    data=(slice(0, 2, None), slice(0, 2, None))),
                              LinOp(type='neg', size=(2, 2),
                                    args=[LinOp(type='dense_const', size=(2, 2), args=[], data=matrix([[ 1.,  0.],
                                                                                                       [ 0.,  1.]]))],
                                    data=None)],
                        data=None), constr_id=3, size=(2, 2)),
 
 LinEqConstr(expr=LinOp(type='sum', size=(2, 1), args=[LinOp(type='index', size=(2, 1), args=[LinOp(type='variable', size=(3, 3), args=[], data=1)], data=(slice(0, 2, None), slice(2, 3, None))), LinOp(type='neg', size=(2, 1), args=[LinOp(type='variable', size=(2, 1), args=[], data=0)], data=None)], data=None), constr_id=4, size=(2, 1)),

 LinEqConstr(expr=LinOp(type='sum', size=(1, 1), args=[LinOp(type='index', size=(1, 1), args=[LinOp(type='variable', size=(3, 3), args=[], data=1)], data=(slice(2, 3, None), slice(2, 3, None))), LinOp(type='neg', size=(1, 1), args=[LinOp(type='variable', size=(1, 1), args=[], data=2)], data=None)], data=None), constr_id=5, size=(1, 1)),
 LinEqConstr(expr=LinOp(type='sum', size=(3, 1), args=[LinOp(type='upper_tri', size=(3, 1), args=[LinOp(type='variable', size=(3, 3), args=[], data=1)], data=None), LinOp(type='neg', size=(3, 1), args=[LinOp(type='upper_tri', size=(3, 1), args=[LinOp(type='transpose', size=(3, 3), args=[LinOp(type='variable', size=(3, 3), args=[], data=1)], data=None)], data=None)], data=None)], data=None), constr_id=7, size=(3, 1)), LinLeqConstr(expr=LinOp(type='neg', size=(6, 1), args=[LinOp(type='mul', size=(6, 1), args=[LinOp(type='reshape', size=(9, 1), args=[LinOp(type='variable', size=(3, 3), args=[], data=1)], data=None)], data=LinOp(type='sparse_const', size=(6, 9), args=[], data=<6x9 sparse matrix of type '<type 'numpy.float64'>'
	with 6 stored elements in Compressed Sparse Column format>))], data=None), constr_id=6, size=(6, 1))]
