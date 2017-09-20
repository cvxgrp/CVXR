library(cvxr)
library(EBImage)
## Load the images.
lena_orig <- channel(readImage("data/lena512.png"), "gray")
lena_corrupt <- channel(readImage("data/lena512_corrupted.png"), "gray")

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
