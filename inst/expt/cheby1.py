import cvxopt
import numpy as np
from pylab import *
import math
import pdb

# from cvxpy import numpy as my_numpy

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
a1 = cvxopt.matrix([2,1], (2,1))
a2 = cvxopt.matrix([2,-1], (2,1))
a3 = cvxopt.matrix([-1,2], (2,1))
a4 = cvxopt.matrix([-1,-2], (2,1))

b = cvxopt.matrix(1, (4,1))


constraints = [ a1.T*center + np.linalg.norm(a1, 2)*radius <= b[0],
                a2.T*center + np.linalg.norm(a2, 2)*radius <= b[1],
                a3.T*center + np.linalg.norm(a3, 2)*radius <= b[2] ]

# objective
objective = Maximize(radius)

p = Problem(objective, constraints)
# The optimal objective is returned by p.solve().
result = p.solve(verbose=True)
# The optimal value
print radius.value
print center.value
