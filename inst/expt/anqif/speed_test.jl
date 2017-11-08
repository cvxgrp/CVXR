using Convex

A = readdlm("a.txt")
b = readdlm("b.txt")
(n, m) = size(A)

x = Variable(m)
problem = minimize(sumsquares(A * x - b))

tic()
solve!(problem)
toc()

problem.status
problem.optval
