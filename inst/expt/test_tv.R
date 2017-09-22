library(cvxr)
if(!("EBImage" %in% installed.packages())) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("EBImage")
}
library(EBImage)
library(imager)


# Load images
orig_img <- load.image("../data/lena512.png")
corr_img <- load.image("..data/lena512_corrupted.png")
orig_img <- flip(resize(orig_img, 64))
corr_img <- flip(resize(corr_img, 64))
colorMode(orig_img) <- "Grayscale"
colorMode(corr_img) <- "Grayscale"
image(orig_img, 1)
image(corr_img, 1)

# Convert to 2-D arrays
Uorig <- imageData(orig_img)[,,1]
Ucorr <- imageData(corr_img)[,,1]
rows <- nrow(Uorig)
cols <- ncol(Uorig)

# Known = 1 if pixel known, 0 if pixel corrupted
Known <- matrix(0, nrow = rows, ncol = cols)
Known[Uorig == Ucorr] <- 1
# image(Uorig, axes = FALSE, col = grey(seq(0, 1, length = 256)))
# image(Ucorr, axes = FALSE, col = grey(seq(0, 1, length = 256)))

U <- Variable(rows, cols)
obj <- Minimize(TotalVariation(U))
constraints <- list(MulElemwise(Known, U) == MulElemwise(Known, Ucorr))
prob <- Problem(obj, constraints)
result <- solve(prob)

Urec <- result$getValue(U)
# image(Urec, axes = FALSE, col = grey(seq(0, 1, length = 256)))
rec_img <- Image(Urec, dim(Urec), "Grayscale")
image(rec_img, 1)
