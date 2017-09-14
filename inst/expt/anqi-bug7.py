from cvxpy import *
import canonInterface
import numpy as np

x = Variable(3)
P = np.eye(3)
obj = matrix_frac(x, P)
constraint = [ x == np.array([1, 5, 3]) ]
prob = Problem(Minimize(obj), constraint)
data = prob.get_problem_data(SCS)
prob.solve(SCS, verbose = True)

##result = prob.solve(SCS, verbose = True)
