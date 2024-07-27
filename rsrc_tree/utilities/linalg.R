## CVXPY SOURCE: cvxpy/utilities/linalg.py

############################
#                          #
# Linear algebra utilities #
#                          #
############################
orth <- function(V, tol = 1e-12, ...) {
  # Return a matrix whose columns are an orthonormal basis for range(V).
  res <- qr(V, tol = tol, ...)
  QR <- res$qr   # Q in upper triangle, R in lower triangle.
  rank <- res$rank
  Q_full <- QR
  Q_full[lower.tri(QR)] <- 0
  Q <- matrix(Q_full[,1:rank], nrow = nrow(V), ncol = rank)   # Ensure 2-dimensional.
  return(Q)
}

onb_for_orthogonal_complement <- function(V) {
  # Let U = the orthogonal complement of range(V).
  # This function returns an array Q whose columns are an orthonormal basis for U.
  # It requires that dim(U) > 0.
  n <- nrow(V)
  Q1 <- orth(V)
  rank <- ncol(Q1)

  if(n <= rank)
    stop("Must have n > rank")

  if(is.complex(V))
    P <- diag(n) - Q1 %*% t(Conj(Q1))
  else
    P <- diag(n) - Q1 %*% t(Q1)
  Q2 <- orth(P)
  return(Q2)
}

is_psd_within_tol <- function(A, tol) {
  # Return TRUE if we can certify that A is PSD (up to tolerance "tol").
  #
  # First we check if A is PSD according to the Gershgorin Circle Theorem.
  #
  # If Gershgorin is inconclusive, then we use an iterative method to estimate
  # extremal eigenvalues of certain shifted versions of A. The shifts are chosen
  # so that the signs of those eigenvalues tell us the signs of the eigenvalues of A.
  #
  # If there are numerical issues then it's possible that this function returns
  # FALSE even when A is PSD. If you know that you're in that situation, then
  # you should replace A by PSDWrap(A).
  #
  # Parameters
  # -----------
  # A = Symmetric (or Hermitian) dense or sparse matrix.
  # tol = Nonnegative floating point variable. Something very small, like 1e-10.

  if(gershgorin_psd_check(A, tol))
    return(TRUE)

  # Returns the eigenvalue w[i] of A where 1/(w[i] - sigma) is minimized.
  # If A - sigma*I is PSD, then w[i] should be equal to the largest eigenvalue of A.
  # If A - sigma*I is not PSD, then w[i] should be the largest eigenvalue of A where w[i] - sigma < 0.
  # We should only call this function with sigma < 0. In this case, if A - sigma*I is not PSD,
  # then A is not PSD, and w[i] < -abs(sigma) is a negative eigenvalue of A.
  # If A - sigma*I is PSD, then we obviously have that the smallest eigenvalue of A is >= sigma.
  SA_eigsh <- function(sigma) {
    require(RSpectra)
    res <- RSpectra::eigs_sym(A, k = 1, which = "SA", sigma = sigma, opts = list(retvec = FALSE))
    return(res$values)
  }

  ev <- NA_real_
  tryCatch({
    ev <- SA_eigsh(-tol)   # Might return NA, or raise an exception.
  }, finally = {
    if(all(is.na(ev))) {
      # Will be NA if A has an eigenvalue which is exactly -tol.
      # (We might also hit this code block for other reasons).
      temp <- tol - .Machine$double.eps
      ev <- SA_eigsh(-temp)
    }
  })
  return(all(ev >= -tol))
}

gershgorin_psd_check <- function(A, tol) {
  # Use the Gershgorin Circle Theorem
  #
  # https://en.wikipedia.org/wiki/Gershgorin_circle_theorem
  #
  # As a sufficient condition for A being PSD with tolerance "tol".
  #
  # The computational complexity of this function is O(nnz(A)).
  #
  # Parameters
  # -----------
  # A = Symmetric (or Hermitian) dense or sparse matrix.
  # tol = Nonnegative floating point variable. Something very small, like 1e-10.

  if(is(A, "sparseMatrix")) {
    d <- diag(A)
    if(any(d < -tol))
      return(FALSE)
    A_shift <- A - sparseMatrix(i = 1:length(d), j = 1:length(d), x = d)
    radii <- apply(A_shift, MARGIN = 2, sum)
    return(all(diag - radii >= -tol))
  } else if(is.numeric(A)) {
    d <- diag(A)
    if(any(diag < -tol))
      return(FALSE)
    A_shift <- A - diag(d)
    A_shift <- abs(A_shift)
    radii <- apply(A_shift, MARGIN = 2, sum)
    return(all(diag - radii >= -tol))
  } else
    stop("A must be a sparse or dense numeric matrix")
}
