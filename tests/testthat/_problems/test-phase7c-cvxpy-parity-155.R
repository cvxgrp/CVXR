# Extracted from test-phase7c-cvxpy-parity.R:155

# test -------------------------------------------------------------------------
skip_if_not_installed("clarabel")
a <- Variable()
b <- Variable()
dom <- domain(kl_div(a, b))
