from cvxpy import *
import numpy as np
import time

m = 500
n = 1000
A = np.loadtxt('a.txt')
b = np.loadtxt('b.txt')

x = Variable(m)
# objective = Minimize(sum_squares(A * x - b))
objective = Minimize(norm(A * x - b))
prob = Problem(objective)

start = time.time()
prob.solve(verbose = True)
print("--- %s seconds ---" % (time.time() - start))

