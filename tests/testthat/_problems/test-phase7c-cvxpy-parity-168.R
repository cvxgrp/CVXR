# Extracted from test-phase7c-cvxpy-parity.R:168

# test -------------------------------------------------------------------------
skip_if_not_installed("clarabel")
a <- Variable()
b <- Variable()
dom <- domain(rel_entr(a, b))
