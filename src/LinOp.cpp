#include "LinOp.hpp"

//' Create a new LinOp object.
//'
//' @return an external ptr (Rcpp::XPtr) to a LinOp object instance.
// [[Rcpp::export(.LinOp__new)]]
SEXP LinOp__new() {
  // create a pointer to an Uniform object and wrap it
  // as an external pointer
  Rcpp::XPtr<LinOp> ptr( new LinOp( ), true );

  // return the external pointer to the R side
  return ptr;
}

//' Get the \code{sparse} flag field for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @return TRUE or FALSE
// [[Rcpp::export(.LinOp__get_sparse)]]
bool LinOp__get_sparse(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  return ptr->sparse;
}

//' Set the flag \code{sparse} of the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @param sparseSEXP an R boolean
// [[Rcpp::export(.LinOp__set_sparse)]]
void LinOp__set_sparse(SEXP xp, SEXP sparseSEXP) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  // Set the result
  ptr->sparse = Rcpp::as<bool>(sparseSEXP);
}


//' Get the field named \code{sparse_data} from the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @return a \link[Matrix]{dgCMatrix-class} object
// [[Rcpp::export(.LinOp__get_sparse_data)]]
Eigen::SparseMatrix<double> LinOp__get_sparse_data(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  return ptr->sparse_data;
}

//' Set the field named \code{sparse_data} of the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @param sparseMat a \link[Matrix]{dgCMatrix-class} object
// [[Rcpp::export(.LinOp__set_sparse_data)]]
void LinOp__set_sparse_data(SEXP xp, SEXP sparseMat) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  // Map to Eigen::SparseMatrix<double> which is Matrix courtesy of Utils.h typedef
  ptr->sparse_data = Rcpp::as<Eigen::SparseMatrix<double> >(sparseMat);
  ptr->sparse = true;

}

//' Get the field \code{dense_data} for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @return a MatrixXd object
// [[Rcpp::export(.LinOp__get_dense_data)]]
Eigen::MatrixXd LinOp__get_dense_data(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  return ptr->dense_data;
}

//' Set the field \code{dense_data} of the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @param denseMat a standard matrix object in R
// [[Rcpp::export(.LinOp__set_dense_data)]]
void LinOp__set_dense_data(SEXP xp, SEXP denseMat) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  // Set the result
  ptr->dense_data = Rcpp::as<Eigen::MatrixXd >(denseMat);
}

//' Get the field \code{size} for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @return an integer vector
// [[Rcpp::export(.LinOp__get_size)]]
std::vector<int>  LinOp__get_size(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  return ptr->size;
}

//' Set the field \code{size} of the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @param value an integer vector object in R
// [[Rcpp::export(.LinOp__set_size)]]
void LinOp__set_size(SEXP xp, Rcpp::IntegerVector value) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  // Set the result
  ptr->size.clear();
  for (int i = 0; i < value.size(); ++i) {
    ptr->size.push_back(value[i]);
  }

}

//' Perform a push back operation on the \code{args} field of LinOp
//'
//' @param xp the LinOp Object XPtr
//' @param yp the LinOp Object XPtr to push
// [[Rcpp::export(.LinOp__args_push_back)]]
void LinOp__args_push_back(SEXP xp, SEXP yp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptrX(xp);
  Rcpp::XPtr<LinOp> ptrY(yp);

  (ptrX->args).push_back(ptrY);
}

//' Perform a push back operation on the \code{size} field of LinOp
//'
//' @param xp the LinOp Object XPtr
//' @param intVal the integer value to push back
// [[Rcpp::export(.LinOp__size_push_back)]]
void LinOp__size_push_back(SEXP xp, int intVal) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  (ptr->size).push_back(intVal);
}

//' Set the field named \code{type} for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @param typeValue an integer value (1-based)
// [[Rcpp::export(.LinOp__set_type)]]
void LinOp__set_type(SEXP xp, int typeValue) {
  OperatorType oType = to_optype(typeValue);
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  ptr->type = oType;
#ifdef _R_DEBUG_
    Rcpp::Rcout << "LinOp Id " << ptr->id << " type now " << ptr->type << std::endl;
#endif
}

//' Get the field named \code{type} for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @return an integer value for type
// [[Rcpp::export(.LinOp__get_type)]]
int LinOp__get_type(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  int result = from_optype(ptr->type);
  return(result);
}

//' Perform a push back operation on the \code{slice} field of LinOp
//'
//' @param xp the LinOp Object XPtr
//' @param intVec an integer vector to push back
// [[Rcpp::export(.LinOp__slice_push_back)]]
void LinOp__slice_push_back(SEXP xp, std::vector<int> intVec) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  (ptr->slice).push_back(intVec);
}

//' Get the slice field of the LinOp Object
//'
//' @param xp the LinOp Object XPtr
//' @return the value of the slice field of the LinOp Object
// [[Rcpp::export(.LinOp__get_slice)]]
std::vector<std::vector<int> >  LinOp__get_slice(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  // Get the result
  return ptr->slice;
}

//' Set the slice field of the LinOp Object
//'
//' @param xp the LinOp Object XPtr
//' @param value a list of integer vectors, e.g. \code{list(1:10, 2L, 11:15)}
//' @return the value of the slice field of the LinOp Object
// [[Rcpp::export(.LinOp__set_slice)]]
void LinOp__set_slice(SEXP xp, std::vector<std::vector<int> > value) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  // set the result
  ptr->slice = value;
}

#ifdef _R_DEBUG
//' Get the id field of the LinOp Object
//'
//' @param xp the LinOp Object XPtr
//' @return the value of the id field of the LinOp Object
// [[Rcpp::export(.LinOp__get_id)]]
std::string LinOp__get_id(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  // Get the result
  return ptr->id;
}
#endif






