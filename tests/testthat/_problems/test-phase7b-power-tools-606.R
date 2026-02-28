# Extracted from test-phase7b-power-tools.R:606

# test -------------------------------------------------------------------------
t_var <- Variable(c(1L, 1L))
x_list <- lapply(1:4, function(i) Variable(c(1L, 1L)))
p <- as.bigq(c(1L, 1L, 1L, 1L), c(4L, 4L, 4L, 4L))
