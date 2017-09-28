# Problem data
m <- 51
L <- 2
h <- L/(m-1)

# Form objective
x <- Variable(m)
y <- Variable(m)
objective <- Minimize(sum(y))

# Form constraints
constraints <- list(x[1] == 0, y[1] == 1, x[m] == 1, y[m] == 1,
                    diff(x)^2 + diff(y)^2 <= h^2)

# Solve catenary problem
prob <- Problem(objective, constraints)
system.time(result <- solve(prob))

# Plot and compare with ideal catenary
xs <- result$getValue(x)
ys <- result$getValue(y)
plot(c(0,1), c(0,1), type = "n", xlab = "x", ylab = "y")
lines(xs, ys, col = "blue", lwd = 2)

points(c(0,1), c(1,1))
curve(0.22964*cosh((x - 0.5)/0.22964) - 0.02603, 0, 1, col = "red", add = TRUE)
grid()
