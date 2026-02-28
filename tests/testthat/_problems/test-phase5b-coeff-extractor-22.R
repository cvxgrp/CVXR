# Extracted from test-phase5b-coeff-extractor.R:22

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
p <- Problem(Minimize(x + 3))
inv_data <- InverseData(p)
