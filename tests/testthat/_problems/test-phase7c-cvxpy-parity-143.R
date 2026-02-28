# Extracted from test-phase7c-cvxpy-parity.R:143

# test -------------------------------------------------------------------------
skip_if_not_installed("clarabel")
a <- Variable()
dom <- domain(log1p_atom(a))
