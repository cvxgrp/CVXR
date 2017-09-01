from cvxpy import *
import numpy as np
import canonInterface
import pandas
#np.random.seed(1)

# Initialize data.
MAX_ITERS = 10
rho = 1.0
n = 20
m = 10
#A = np.random.randn(m,n)
#b = np.random.randn(m,1)
A = pandas.read_csv("a.csv", header=None).values
b = pandas.read_csv("b.csv", header=None).values
# Initialize problem.
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


soln = Problem(Minimize(aug_lagr)).solve()

for t in range(MAX_ITERS):
    Problem(Minimize(aug_lagr)).solve()
    y.value += rho*resid.value
    
print "Optimal value from method of multipliers", f.value
