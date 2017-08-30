#! /usr/bin/env python

from cvxpy import *
import canonInterface
import numpy

n = 2
A = numpy.array([[1, 0],
                 [0, 1]])
x = Variable(n)
objective = Minimize(sum(A * x))
constraint = [x >= 1]
prob = Problem(objective, constraint)

result = prob.solve(verbose = True)
