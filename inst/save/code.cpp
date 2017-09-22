
// [[Rcpp::export]]
Rcpp::IntegerMatrix IndexGetMat(Rcpp::IntegerVector expr_size,
				SEXP row,
				Rcpp::IntegerVector col) {
  int expr_size_0 = expr_size[0], expr_size_1 = expr_size[1];
  int expr_prod = expr_size_0 * expr_size_1;
  int *v = new int[expr_prod];
  
  for (int i = 0; i < expr_prod; ++i) {
    v[i] = i;
  }
  
  Rcpp::IntegerVector select_mat;
  
  if (Rf_isMatrix(row) && Rf_NilValue(col)) {
    // Do matrix stuff since row is a matrix
    longint row_size = row
    select_mat <- Rcpp::IntegerVector(
    
  } else if (Rf_NilValue(row) && Rf_NilValue(col)) {
    //
  } else if (!Rf_NilValue(row) && Rf_NilValue(col)) {
    //
  } else {
    
  }

    delete [] v;
    
  v.attr("dim") = Dimension(2, 2);
    
    NumericVector rcpp_matrix(){
    // Creating a vector object
    NumericVector v = {1,2,3,4};

    // Set number of rows and columns to attribute dim
    v.attr("dim") = Dimension(2, 2);

    // Returns a vector with attribute dim to R
    return v;
}
  
}
				
