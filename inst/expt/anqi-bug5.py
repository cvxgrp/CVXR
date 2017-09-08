from cvxpy import *
import canonInterface
import numpy as np

x = Variable(2)
P = np.eye(2)
obj = matrix_frac(x, P)
prob = Problem(Minimize(obj))
data = prob.get_problem_data(SCS)

##result = prob.solve(SCS, verbose = True)

