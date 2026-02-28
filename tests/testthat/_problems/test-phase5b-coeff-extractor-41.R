# Extracted from test-phase5b-coeff-extractor.R:41

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
p <- Problem(Minimize(2 * x))
inv_data <- InverseData(p)
