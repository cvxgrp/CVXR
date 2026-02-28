# Extracted from test-phase5b-coeff-extractor.R:12

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(3)
p <- Problem(Minimize(sum_entries(x)))
inv_data <- InverseData(p)
