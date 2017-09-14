from cvxpy import *
import numpy as np
import canonInterface

rho = 1.0
n = 3
m = 2
A = np.array([[-0.5604756, 1.55870831, 0.1292877],
              [-0.2301775, 0.07050839, 1.7150650]])
b = np.array([0.4609162, -1.2650612])

x = Variable(n)
f = norm(x, 1)

# Solve with CVXPY.
p = Problem(Minimize(f), [A*x == b])
p.solve()
print "Optimal value from CVXPY", f.value

# Solve with method of multipliers.
resid = A*x - b
y = Parameter(m); y.value = np.zeros(m)
aug_lagr = f + y.T*resid + (rho/2)*sum_squares(resid)
p = Problem(Minimize(aug_lagr))
p.solve()

# b canonInterface.set_matrix_data
#data = p.get_problem_data(ECOS)
# print "Optimal value from method of multipliers", f.value
