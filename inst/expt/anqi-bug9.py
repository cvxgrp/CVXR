from cvxpy import *
import canonInterface

x_bool = Bool()
obj = Minimize(square(x_bool - 0.2))
p = Problem(obj, list())
data = p.get_problem_data(ECOS_BB)

soln = p.solve(solver = ECOS_BB)


