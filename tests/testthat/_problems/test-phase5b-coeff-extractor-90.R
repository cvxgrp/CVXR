# Extracted from test-phase5b-coeff-extractor.R:90

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(1)
p <- Problem(Minimize(x))
inv_data <- InverseData(p)
