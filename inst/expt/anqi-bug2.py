#! /usr/bin/env python

from cvxpy import *
import canonInterface
import numpy
import scipy as sp

x = Variable(1)
obj = Maximize(log(x))
prob = Problem(obj)
sol = prob.solve(verbose = True)
data = prob.get_problem_data(ECOS)

d = numpy.matrix([[ 0., -1.],
                  [-1.,  0.],
                  [ 0.,  0.]])

c = numpy.array([ 0., -1.])
G = sp.sparse.csc_matrix(d)
h = numpy.array([-0., -0.,  1.])
dims = {'ep': 1, 'f': 0, 'l': 0, 'e': 1, 'q': [], 's': []}


d1 = numpy.matrix([[ 0., -1.],
                   [ 0.,  0.],
                   [-1.,  0.]])

c1 = numpy.array([ 0., -1.])
G1 = sp.sparse.csc_matrix(d1)
h1 = numpy.array([-0., 1,  0])
dims1 = {'ep': 1, 'f': 0, 'l': 0, 'e': 1, 'q': [], 's': []}

ecos.solve(c1,G1,h1,dims1, verbose = True)

