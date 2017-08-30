In [13]: prob.get_problem_data(ECOS)
Out[13]: 
{'A': <0x2 sparse matrix of type '<type 'numpy.float64'>'
 	with 0 stored elements in Compressed Sparse Column format>,
 'F': None,
 'G': <3x2 sparse matrix of type '<type 'numpy.float64'>'
 	with 2 stored elements in Compressed Sparse Column format>,
 'b': array([], dtype=float64),
 'bool_vars_idx': [],
 'c': array([ 0., -1.]),
 'dims': {'ep': 1, 'f': 0, 'l': 0, 'q': [], 's': []},
 'h': array([-0., -0.,  1.]),
 'int_vars_idx': [],
 'offset': 0.0}

In [14]: data = prob.get_problem_data(ECOS)

In [15]: data['A'].todense
Out[15]: 
<bound method csc_matrix.todense of <0x2 sparse matrix of type '<type 'numpy.float64'>'
	with 0 stored elements in Compressed Sparse Column format>>

In [16]: data['A'].todense()
Out[16]: matrix([], shape=(0, 2), dtype=float64)

In [17]: data['G'].todense()
Out[17]: 
matrix([[ 0., -1.],
        [-1.,  0.],
        [ 0.,  0.]])

In [18]: 
