import numpy as np
from cvxpy import *

# Taken from CVX website http://cvxr.com/cvx/examples/
# Example: Compute and display the Chebyshev center of a 2D polyhedron
# Ported from cvx matlab to cvxpy by Misrab Faizullah-Khan
# Original comments below

# Boyd & Vandenberghe, "Convex Optimization"
# Joelle Skaf - 08/16/05
# (a figure is generated)
#
# The goal is to find the largest Euclidean ball (i.e. its center and
# radius) that lies in a polyhedron described by linear inequalites in this
# fashion: P = { x : a_i'*x <= b_i, i=1,...,m } where x is in R^2

# Create the problem

# variables
radius = Variable(1)
center = Variable(2)

# constraints
A = np.array([[2, 2, -1, -1],
              [1, -1, 2, -2]])

a1 = np.array([[2, 1]])
a2 = np.array([[2,-1]])
a3 = np.array([[-1, 2]])
a4 = np.array([[-1,-2]])

b = np.array([1, 1, 1, 1])

constraints = [ a1 * center + norm(a1, 2) * radius <= b[0],
                a2 * center + norm(a2, 2) * radius <= b[1],
                a3 * center + norm(a3, 2) * radius <= b[2],
                a4 * center + norm(a4, 2) * radius <= b[3] ]

# objective
objective = Maximize(radius)

p = Problem(objective, constraints)
# The optimal objective is returned by p.solve().
result = p.solve(verbose=True)
# The optimal value
print radius.value
print center.value
