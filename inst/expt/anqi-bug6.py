from cvxpy import *
import canonInterface
import numpy as np
import ecos

x = np.array([[1],
              [-1]])
p = Problem(Minimize(max_elemwise(x.T, 2, 2 + x.T)[1]))

from cvxpy.problems.solvers.ecos_intf import ECOS
# b ECOS.solve

p.solve()

{'A': <0x0 sparse matrix of type '<type 'numpy.float64'>'
	with 0 stored elements in Compressed Sparse Column format>, 'c': array([], dtype=float64), 'b': array([], dtype=float64), 'G': <0x0 sparse matrix of type '<type 'numpy.float64'>'
	with 0 stored elements in Compressed Sparse Column format>, 'F': None, 'bool_vars_idx': [], 'h': array([], dtype=float64), 'dims': {'q': [], 's': [], 'f': 0, 'l': 0, 'ep': 0}, 'offset': 2.0, 'int_vars_idx': []}

