# Extracted from test-phase5b-coeff-extractor.R:73

# prequel ----------------------------------------------------------------------
library(CVXR)

# test -------------------------------------------------------------------------
x <- Variable(2)
p <- Problem(Minimize(sum_entries(x)), list(x >= 0))
inv_data <- InverseData(p)
