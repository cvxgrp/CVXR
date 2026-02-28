## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
set.seed(123)
n <- 100
p <- 10
beta <- -4:5

X <- matrix(rnorm(n * p), nrow = n)
Y <- X %*% beta + rnorm(n)

## -----------------------------------------------------------------------------
ls.model <- lm(Y ~ 0 + X)

## -----------------------------------------------------------------------------
library(CVXR)

betaHat <- Variable(p)
objective <- Minimize(sum((Y - X %*% betaHat)^2))
problem <- Problem(objective)
result <- psolve(problem)

## -----------------------------------------------------------------------------
cat("Optimal value:", result, "\n")
cbind(CVXR = round(value(betaHat), 3),
      lm   = round(coef(ls.model), 3))

## -----------------------------------------------------------------------------
problem <- Problem(objective, constraints = list(betaHat >= 0))
result <- psolve(problem)
round(value(betaHat), 3)

## -----------------------------------------------------------------------------
A <- matrix(c(0, 1, 1, rep(0, 7)), nrow = 1)
B <- diag(c(1, 0, 0, rep(1, 7)))

constraint1 <- A %*% betaHat <= 0
constraint2 <- B %*% betaHat >= 0

problem <- Problem(objective, constraints = list(constraint1, constraint2))
result <- psolve(problem, verbose = TRUE) ## verbose = TRUE for details
round(value(betaHat), 3)

## -----------------------------------------------------------------------------
installed_solvers()

## ----eval = FALSE-------------------------------------------------------------
# psolve(problem, solver = "CLARABEL")

## -----------------------------------------------------------------------------
sessionInfo()

