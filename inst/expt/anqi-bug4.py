from cvxpy import *
import canonInterface

A = Variable(2,2)
obj = Maximize(log_det(A))
p = Problem(obj)

result = p.solve(SCS, verbose = True)

