# Load images
data(lena)
data(lena_corr)
rows <- nrow(lena)
cols <- ncol(lena)

# Known = 1 if pixel known, 0 if pixel corrupted
Known <- matrix(0, nrow = rows, ncol = cols)
Known[lena == lena_corr] <- 1
image(lena, axes = FALSE, col = grey(seq(0, 1, length = 256)))
image(lena_corr, axes = FALSE, col = grey(seq(0, 1, length = 256)))

U <- Variable(rows, cols)
obj <- Minimize(TotalVariation(U))
constraints <- list(Known * U == Known * lena_corr)
prob <- Problem(obj, constraints)
result <- solve(prob)
# TODO: SOC_AXIS LinOp is not implemented?

# TODO: More user-friendly functions to retrieve results
result$optimal_value
lena_recon <- result$primal_values[[as.character(U@id)]]
image(lena_recon, axes = FALSE, col = grey(seq(0, 1, length = 256)))