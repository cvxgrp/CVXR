x_bool <- Bool()
y_int <- Int()
A_bool <- Bool(3,2)
B_int <- Int(2,3)

test_that("Test that MIP problems are deterministic", {
  data_recs <- list()
  result_recs <- list()
  for(i in 1:5) {
    obj <- Minimize(Square(y_int - 0.2))
    p <- Problem(obj, list(A_bool == 0, x_bool == B_int))
    # data_recs <- c(data_recs, get_problem_data(p, "ECOS_BB"))
  }
  
  # Check that problem data and result is always the same
  for(i in 1:5) {
    for(key in c("c", "A", "b", "G", "h", "bool_vars_idx", "int_vars_idx")) {
      lh_item <- data_recs[1][key]
      rh_item <- data_recs[i][key]
      if(key %in% c("A", "G")) {
        lh_item <- as.matrix(lh_item)
        rh_item <- as.matrix(rh_item)
      }
      expect_equal(lh_item, rh_item, tolerance = TOL)
    }
  }
})