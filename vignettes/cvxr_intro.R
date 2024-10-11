## --------------------------------------------------------------------------------
set.seed(123)
n <- 100
p <- 10
beta <- -4:5   # beta is just -4 through 5.
X <- matrix(rnorm(n * p), nrow=n)
colnames(X) <- paste0("beta_", beta)
Y <- X %*% beta + rnorm(n)
ls.model <- lm(Y ~ 0 + X)   # There is no intercept in our model above
m <- data.frame(ls.est = coef(ls.model))
rownames(m) <- paste0("$\\beta_{", 1:p, "}$")
knitr::kable(m)

suppressWarnings(library(CVXR, warn.conflicts=FALSE))

betaHat <- Variable(p)

objective <- Minimize(sum((Y - X %*% betaHat)^2))

problem <- Problem(objective)

result <- solve(problem)


## ----echo = FALSE----------------------------------------------------------------
solution <- result$getValue(betaHat)
cat(sprintf("Objective value: %f\n", result$value))


## --------------------------------------------------------------------------------
m <- cbind(coef(ls.model), result$getValue(betaHat))
colnames(m) <- c("lm est.", "CVXR est.")
rownames(m) <- paste0("$\\beta_{", 1:p, "}$")
knitr::kable(m)


## --------------------------------------------------------------------------------
problem <- Problem(objective, constraints = list(betaHat >= 0))
result <- solve(problem)
m <- data.frame(CVXR.est = result$getValue(betaHat))
rownames(m) <- paste0("$\\beta_{", 1:p, "}$")
knitr::kable(m)


## --------------------------------------------------------------------------------
if (requireNamespace("nnls", quietly = TRUE)) {
    nnls.fit <- nnls::nnls(X, Y)$x
} else {
    nnls.fit <- rep(NA, p)
}


## --------------------------------------------------------------------------------
m <- cbind(result$getValue(betaHat), nnls.fit)
colnames(m) <- c("CVXR est.", "nnls est.")
rownames(m) <- paste0("$\\beta_{", 1:p, "}$")
knitr::kable(m)


## --------------------------------------------------------------------------------
A <- matrix(c(0, 1, 1, rep(0, 7)), nrow = 1)
colnames(A) <- paste0("$\\beta_{", 1:p, "}$")
knitr::kable(A)


## --------------------------------------------------------------------------------
constraint1 <- A %*% betaHat <= 0


## ----eval = FALSE----------------------------------------------------------------
## constraint1 <- betaHat[2] + betaHat[3] <= 0


## --------------------------------------------------------------------------------
B <- diag(c(1, 0, 0, rep(1, 7)))
colnames(B) <- rownames(B) <- paste0("$\\beta_{", 1:p, "}$")
    knitr::kable(B)


## --------------------------------------------------------------------------------
constraint2 <- B %*% betaHat >= 0


## --------------------------------------------------------------------------------
problem <- Problem(objective, constraints = list(constraint1, constraint2))
result <- solve(problem)


## --------------------------------------------------------------------------------
m <- data.frame(CVXR.soln = result$getValue(betaHat))
rownames(m) <- paste0("$\\beta_{", 1:p, "}$")
knitr::kable(m)

