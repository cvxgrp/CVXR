from cvxpy import *
import canonInterface
#from cvxpy.atoms.affine import conv

n = 3
x = Variable(n)
f = [1, 2, 3]
g = [0, 1, 0.5]
f_conv_g = [0., 1., 2.5,  4., 1.5]
expr = conv(f, x)

# Matrix stuffing.
t = Variable()
prob = Problem(Minimize(norm(expr, 1)),
               [x == g])
result = prob.solve()
# b conv.value
w = expr.value

##self.assertAlmostEqual(result, sum(f_conv_g))
##self.assertItemsAlmostEqual(expr.value, f_conv_g)
