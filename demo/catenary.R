# Problem data
m <- 101
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
result <- solve(prob)

# Plot and compare with ideal catenary
xs <- result$getValue(x)
ys <- result$getValue(y)
plot(c(0,1), c(0,1), type = "n", xlab = "x", ylab = "y")
lines(xs, ys, col = "blue", lwd = 2)

points(c(0,1), c(1,1))
curve(0.22964*cosh((x - 0.5)/0.22964) - 0.02603, 0, 1, col = "red", add = TRUE)
grid()

# Lower right endpoint and add staircase structure
ground <- sapply(seq(0, 1, length.out = m), function(x) {
  if(x < 0.2)
    return(0.6)
  else if(x >= 0.2 && x < 0.4)
    return(0.4)
  else if(x >= 0.4 && x < 0.6)
    return(0.2)
  else
    return(0)
})

# Solve catenary problem with ground constraint
constraints <- c(constraints, y >= ground)
constraints[[4]] <- (y[m] == 0.5)
prob <- Problem(objective, constraints)
result <- solve(prob)

# Plot catenary against ground
xs <- result$getValue(x)
ys <- result$getValue(y)
plot(c(0, 1), c(1, 0.5), type = "n", xlab = "x", ylab = "y", ylim = c(0, 1))
points(c(0, 1), c(1, 0.5))
lines(xs, ys, col = "blue", lwd = 2)
lines(xs, ground, col = "red")
grid()
