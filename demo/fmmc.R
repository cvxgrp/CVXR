# Boyd, Diaconis, and Xiao. SIAM Rev. 46 (2004) pgs. 667-689 at pg. 672
# Form the complementary graph
antiadjacency <- function(g) {
  n <- max(as.numeric(names(g)))   # Assumes names are integers starting from 1
  a <- lapply(1:n, function(i) c())
  names(a) <- 1:n
  for(x in names(g)) {
    for(y in 1:n) {
      if(!(y %in% g[[x]]))
        a[[x]] <- c(a[[x]], y)
    }
  }
  a
}

# Fastest mixing Markov chain on graph g
FMMC <- function(g, verbose = FALSE) {
  a <- antiadjacency(g)
  n <- length(names(a))
  P <- Variable(n, n)
  o <- rep(1, n)
  objective <- Minimize(norm(P - 1.0/n, "2"))
  constraints <- list(P %*% o == o, t(P) == P, P >= 0)
  for(i in names(a)) {
    for(j in a[[i]]) {  # (i-j) is a not-edge of g!
      idx <- as.numeric(i)
      if(idx != j)
        constraints <- c(constraints, P[idx,j] == 0)
    }
  }
  prob <- Problem(objective, constraints)
  result <- solve(prob)
  if(verbose)
    cat("Status: ", result$status, ", Optimal Value = ", result$value)
  list(status = result$status, value = result$value, P = result$getValue(P))
}

disp_result <- function(states, P, tol = 1e-3) {
  if(require("markovchain")) {
    P[P < tol] <- 0
    P <- P/apply(P, 1, sum)   # Normalize so rows sum to exactly 1
    mc <- new("markovchain", states = states, transitionMatrix = P)
    plot(mc)
  } else {
    rownames(P) <- states
    colnames(P) <- states
    print(P)
  }
}

# SIAM Rev. 46 examples pg. 674: Figure 1 and Table 1
# a) line graph L(4)
g <- list("1" = 2, "2" = c(1,3), "3" = c(2,4), "4" = 3)
result <- FMMC(g, verbose = TRUE)
disp_result(names(g), result$P)

# b) triangle + one edge
g <- list("1" = 2, "2" = c(1,3,4), "3" = c(2,4), "4" = c(2,3))
result <- FMMC(g, verbose = TRUE)
disp_result(names(g), result$P)

# c) bipartite 2 + 3
g <- list("1" = c(2,4,5), "2" = c(1,3), "3" = c(2,4,5), "4" = c(1,3), "5" = c(1,3))
result <- FMMC(g, verbose = TRUE)
disp_result(names(g), result$P)

# d) square + central point
g <- list("1" = c(2,3,5), "2" = c(1,4,5), "3" = c(1,4,5), "4" = c(2,3,5), "5" = c(1,2,3,4,5))
result <- FMMC(g, verbose = TRUE)
disp_result(names(g), result$P)
