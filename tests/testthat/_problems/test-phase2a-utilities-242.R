# Extracted from test-phase2a-utilities.R:242

# test -------------------------------------------------------------------------
caught <- FALSE
tryCatch(
    stop(SolverError("solver failed")),
    SolverError = function(e) {
      caught <<- TRUE
    }
  )
