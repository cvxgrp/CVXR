APPROX_CONES <- list("RelEntrConeQuad" = list("SOC"),
                     "OpRelEntrConeQuad" = list("PSD"))

# Helper function for returning the weights and nodes for an n-point Gauss-Legendre quadrature on [0, 1].
gauss_legendre <- function(n) {
  beta <- 0.5/sqrt(rep(1, n-1) - (2*seq(1, n-1))^(-2))
  
  Tmat <- matrix(0, nrow = n, ncol = n)
  for(i in seq(1, n - 1)) {
    Tmat[i,i+1] <- beta[i]
    Tmat[i+1,i] <- beta[i]
  }
  
  Tmat_eig <- eigen(Tmat, only.values = FALSE)
  D <- Tmat_eig$values
  V <- Tmat_eig$vectors
  
  x <- D
  x <- sort(x)
  i <- order(x)
  w <- 2*(V[1,i])^2
  x <- (x + 1)/2
  w <- w/2
  return(list(w, x))
}
