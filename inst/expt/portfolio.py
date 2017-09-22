import numpy as np
##np.random.seed(123)
n = 10
# mu = np.abs(np.random.randn(n, 1))
# Sigma = np.random.randn(n, n)
# Sigma = Sigma.T.dot(Sigma)

Sigma = np.genfromtxt("Sigma.csv", delimiter=",")
mu = np.genfromtxt("mu.csv", delimiter=",")

# Long only portfolio optimization.
from cvxpy import *
w = Variable(n)
gamma = Parameter(sign='positive')
gamma.value = 1
ret = mu.T*w 
risk = quad_form(w, Sigma)

constraints =  [sum_entries(w) == 1, 
                w >= 0]

prob = Problem(Maximize(ret - gamma*risk),
               constraints)
soln = prob.solve(verbose = True)

constraints[0].dual_value
constraints[1].dual_value
