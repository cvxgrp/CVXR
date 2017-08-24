#! /usr/bin/env python
from cvxpy import *
import canonInterface
import numpy

x  = numpy.random.randn(4, 10)
A = numpy.cov(x)
#S = Semidef(4)
S = A

delta = 0.1
alpha = Variable(4)

e1 = numpy.array([1, 0, 0, 0])
e2 = numpy.array([0, 1, 0, 0])
e3 = numpy.array([0, 0, 1, 0])
e4 = numpy.array([0, 0, 0, 1])

S = numpy.genfromtxt("s.csv", delimiter=",")

objective = Minimize(0.5 * quad_form(alpha, S) + alpha.T * e2 + delta * norm(alpha, 1))
#constraints = [ alpha.T * S * alpha >= 0]
prob = Problem(objective)
prob.solve()
