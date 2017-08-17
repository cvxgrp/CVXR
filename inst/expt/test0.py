#! /usr/bin/env python
from cvxpy import *
import canonInterface
import numpy

# Problem data.
n = 2

A = numpy.array([[1, 0],
                 [0, 1]])
# Construct the problem.
x = Variable(n)
objective = Minimize(sum_entries(A*x))
w = objective.canonicalize()
constraints = [1 <= x]
prob = Problem(objective, constraints)

# Set breakpoint at build_matrix of CVXcanon
##b CVXcanon.build_matrix
## b canonInterface.get_problem_matrix
# The optimal objective is returned by prob.solve().
result = prob.solve()
# The optimal value for x is stored in x.value.
print x.value
# The optimal Lagrange multiplier for a constraint
# is stored in constraint.dual_value.
print constraints[0].dual_value
