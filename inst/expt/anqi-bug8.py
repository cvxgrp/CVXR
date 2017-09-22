from cvxpy import *
import canonInterface

# This throws an error: Subscript out of bounds in saveValuesById
#a = Variable(1)
#obj = Minimize(0*a)
#p = Problem(obj)
#result = p.solve()

# This runs, but throws a warning: A->p (column pointers) not strictly increasing, column 24 empty.
x = Variable(2,2)
obj = Minimize(tv(x))
prob = Problem(obj)
result = prob.solve("SCS")
