library(cvxr)
library(png)

# Load images
orig_img <- readPNG("data/lena512.png")
corr_img <- readPNG("data/lena512_corrupted.png")

# Convert to 2-D arrays
rotate <- function(x) t(apply(x, 2, rev))
Uorig <- rotate(orig_img[,,1])
Ucorr <- rotate(corr_img[,,1])
rows <- nrow(Uorig)
cols <- ncol(Uorig)

# Known = 1 if pixel known, 0 if pixel corrupted
Known <- matrix(0, nrow = rows, ncol = cols)
Known[Uorig == Ucorr] <- 1
image(Uorig, axes = FALSE, col = grey(seq(0, 1, length = 256)))
image(Ucorr, axes = FALSE, col = grey(seq(0, 1, length = 256)))

U <- Variable(rows, cols)
obj <- Minimize(TotalVariation(U))
constraints <- list(MulElemwise(Known, U) == MulElemwise(Known, Ucorr))
prob <- Problem(obj, constraints)
result <- solve(prob)
# TODO: SOC_AXIS LinOp is not implemented?