import cvxpy as cp
import numpy as np
import pickle as pk
with open("xybeta.pickle", "rb") as f:
    l = pk.load(f)
X = l['X']
Y = l['Y'].flatten()
beta = l['beta']

n = X.shape[0]
p = X.shape[1]

betaHat = cp.Variable(p)

cost = cp.sum_squares(Y - X @ betaHat)
problem = cp.Problem(cp.Minimize(cost))

problem.solve(verbose = True)

