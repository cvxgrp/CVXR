#! /usr/bin/env python

from cvxpy import *
import canonInterface
import numpy

x = Variable(1)
obj = Maximize(log(x))
prob = Problem(obj)
sol = prob.solve()
