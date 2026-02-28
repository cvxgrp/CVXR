# Extracted from test-phase2c-operators.R:348

# test -------------------------------------------------------------------------
expect_error(MulExpression(Variable(c(3, 4)), Variable(c(5, 2))),
               "Incompatible")
