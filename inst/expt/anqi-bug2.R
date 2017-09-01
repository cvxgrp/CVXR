library(cvxr)
x <- Variable(1)
obj <- Maximize(Log(x))
prob <- Problem(obj)
data <- get_problem_data(prob, "ECOS")
base::trace("Solver.solve", tracer = browser, exit = browser, signature = c("ECOS"))
solve(prob, verbose = 1)

cat("c is \n")
print(data$c)
cat("G is \n")
print(data$G)
cat("h is \n")
print(data$h)
cat("dims is\n")
print(data$dims)

data$dims$e <- data$dims$ep
ns <- names(data$dims)
data$dims <- lapply(data$dims, as.integer)
names(data$dims) <- ns

ECOSolveR::ECOS_csolve(data$c, data$G, data$h, data$dims, control = ecos.control(verbose=1L))

## Correct answer, i.e. what cvxpy produces for G, and h. G and h are the only difference!
G1 <- Matrix(c(0., -1.,
               -1.,  0.,
               0.,  0), byrow = TRUE, ncol = 2, sparse = TRUE)

h1 = c(0, 0, 1.0)

ECOSolveR::ECOS_csolve(data$c, G1, h1, data$dims, control = ecos.control(verbose=1L))
