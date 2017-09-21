library(cvxr)
library(EBImage)
## Load the images.
lena_orig <- channel(readImage("data/lena512.png"), "gray")
lena_corrupt <- channel(readImage("data/lena512_corrupted.png"), "gray")

##
lena_orig <- resize(lena_orig, 16, 16)
lena_corrupt <- resize(lena_corrupt, 16, 16)
# Convert to arrays.
dims <- dim(lena_orig)

## Known is 1 if the pixel is known,
## 0 if the pixel was corrupted.
known <- (lena_orig@.Data == lena_corrupt@.Data) * 1

U <- Variable(dims[1], dims[2])
obj <- Minimize(TotalVariation(U))
constraints <- [mul_elemwise(Known, U) == mul_elemwise(Known, Ucorr)]
prob <- Problem(obj, constraints)
# Use SCS to solve the problem.
solve(prob, verbose=True, solver="SCS")


library(cvxr)
library(png)

# Load images
orig_img <- readPNG("data/lena512.png")
corr_img <- readPNG("data/lena512_corrupted.png")

# Convert to 2-D arrays
rotate <- function(x) t(apply(x, 2, rev))
Uorig <- rotate(orig_img[,,1])
Ucorr <- rotate(corr_img[,,1])
Uorig <- resize(Uorig, 128, 128)
Ucorr <- resize(Ucorr, 128, 128)

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
result <- solve(prob, solver="SCS", verbose = TRUE)

