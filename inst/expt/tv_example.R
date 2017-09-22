library(cvxr)
## if(!("EBImage" %in% installed.packages())) {
##     source("https://bioconductor.org/biocLite.R")
##     biocLite("EBImage")
## }
## library(EBImage)

                                        # Load images
library(imager)
orig_img <- imager::load.image("data/lena512.png")
corr_img <- imager::load.image("data/lena512_corrupted.png")
orig_img <- imager::grayscale(orig_img)
corr_img <- imager::grayscale(corr_img)

orig_img <- imager::resize(orig_img, 64, 64)
corr_img <- imager::resize(corr_img, 64, 64)

plot(orig_img, axes=FALSE)
plot(corr_img, axes=FALSE)

# Convert to 2-D arrays
Uorig <- orig_img[, , 1, 1]
Ucorr <- corr_img[, , 1, 1]
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
rec_img <- as.cimg(Urec)

plot(rec_img, axes = FALSE)
