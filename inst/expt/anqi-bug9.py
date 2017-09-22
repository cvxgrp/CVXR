from cvxpy import *
import canonInterface

# This throws an error: Subscript out of bounds in saveValuesById
#a = Variable(1)
#obj = Minimize(0*a)
#p = Problem(obj)
#result = p.solve()
x_bool = Bool()
obj = Minimize(square(x_bool - 0.2))
p = Problem(obj, list())
data = p.get_problem_data(ECOS_BB)

#base::trace("Solver.solve", tracer = browser, exit = browser, signature = c("ECOS_BB"))
soln = solve(p, solver = "ECOS_BB")



# This runs, but throws a warning: A->p (column pointers) not strictly increasing, column 24 empty.
x = Variable(2,2)
obj = Minimize(tv(x))
prob = Problem(obj)
result = prob.solve("SCS")
