# Extracted from test-phase2a-utilities.R:224

# test -------------------------------------------------------------------------
caught <- FALSE
tryCatch(
    stop(DCPError("dcp violation")),
    DCPError = function(e) {
      caught <<- TRUE
      expect_equal(e$message, "dcp violation")
    }
  )
