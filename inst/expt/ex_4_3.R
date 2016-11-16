library(cvxr)

# Problem data.


n <- 3

P <- matrix(c(13, 12, -2,
              12, 17, 6,
              -2, 6, 12),
            byrow = TRUE, nrow = n)

q <- t(matrix(c(-22, -14.5, 13), nrow = n))
r <- 1
x_star <- matrix(c(1, 1/2, -1), nrow = n)

# Frame and solve the problem

x <- Variable(n)
objective <- Minimize(  0.5 * quad_form(x, P)  + q.T * x + r )
constraints = list( x >= -1, x <= 1)

base::trace("solve", tracer=browser, exit = browser, signature = c("Problem"))

##debug(solve)

##debug(build_lin_op_tree)

w <- solve(prob)

