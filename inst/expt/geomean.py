from cvxpy import *
import canonInterface

x = Variable(2)
cost = geo_mean(x)
prob = Problem(Maximize(cost), [x <= 1])
result = prob.solve(verbose = True)
