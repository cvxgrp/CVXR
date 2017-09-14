from cvxpy import *
import numpy as np
import canonInterface
from cvxpy.problems.solvers.ecos_intf import ECOS
m = 3
n = 2
# A = np.array([[-0.560475646552213,0.070508391424576],
#               [-0.23017748948328,0.129287735160946],
#               [1.55870831414912,1.71506498688328]])

# b = np.array([[0.460916205989202],
#               [-1.26506123460653],
#               [-0.686852851893526]])
A = np.array([[1, 1],
              [-1, 2],
              [-2, -3]])

b = np.array([[1],
              [2],
              [3]])

## Construct the problem
x = Variable(n)
objective = Minimize(sum_squares(A * x - b))
constraints = [x >= 1, x <= 2]
p = Problem(objective, constraints)

## 
## b ECOS.solve
## The optimal objective is returned by solve(p)
result = p.solve(verbose = True)

