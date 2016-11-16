import cvxopt
import numpy as np
from pylab import *
import math
import pdb

# from cvxpy import numpy as my_numpy

from cvxpy import *

n = 5
L = 2
h = L/(n-1)
h2 = h*h
x = Variable(2*n)

##B = np.diag(np.ones(2*n))a
##B[range(n), range(n)] = 0

objective = Minimize(sum_entries(x[n:]))

constraints = [ x >= 0,
                x[0] == 0,
                x[n-1] == 1,
                x[n] == 1,
                x[2*n-1] == 1 ]

for i in range(n-1):
    A = np.zeros((2, 2*n))
    A[0, i] = -1
    A[0, i+1] = 1
    A[1, n+i] = -1
    A[1, n+i+1] = 1
    constraints += [ norm(A*x) <= h ]

p = Problem(objective, constraints)
# The optimal objective is returned by p.solve().
result = p.solve(verbose=True)

