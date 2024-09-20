#include "LinOp.hpp"
#include "LinOpVector.hpp"
#include "Utils.h"

//' Create a new LinOp object.
//'
//' @return an external ptr (Rcpp::XPtr) to a LinOp object instance.
// [[Rcpp::export(.LinOp__new)]]
Rcpp::XPtr<LinOp> LinOp__new(int typeValue, const std::vector<int> &shape,
		const Rcpp::XPtr<ConstLinOpVector> argsVector) {
  // create a pointer to an LinOp object and wrap it
  // as an external pointer
  Rcpp::XPtr<LinOp> ptr( new LinOp(to_optype(typeValue), shape, argsVector->linop_vec), true );
  
  // return the external pointer to the R side
  return ptr;
}

//' Get the field named \code{type} for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @return an integer value for type
// [[Rcpp::export(.LinOp__get_type)]]
int LinOp__get_type(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  return from_optype(ptr->get_type());
}

//' Return a flag indicating if LinOp is a constant
//'
//' @param xp the LinOp Object XPtr
//' @return an boolean value
// [[Rcpp::export(.LinOp__is_constant)]]
bool LinOp__is_constant(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  return ptr->is_constant();
}

//' Get the \code{dim} field for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @return an integer vector
// [[Rcpp::export(.LinOp__get_dim)]]
std::vector<int>  LinOp__get_dim(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  return ptr->get_shape();
}

//' Get the \code{args} field for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @return a sexp for the args
// [[Rcpp::export(.LinOp__get_args)]]
Rcpp::XPtr<ConstLinOpVector> LinOp__get_args(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> linop_ptr(xp);
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<ConstLinOpVector> ptr(new ConstLinOpVector(), true);
  ptr->linop_vec = linop_ptr->get_args();
  return ptr;
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
  return ptr->get_slice();
}

//' Perform a push back operation on the \code{slice} field of LinOp
//'
//' @param xp the LinOp Object XPtr
//' @param intVec an integer vector to push back
// [[Rcpp::export(.LinOp__push_back_slice_vec)]]
void LinOp__push_back_slice_vec(SEXP xp, std::vector<int> slice_vec) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  ptr->push_back__slice_vec(slice_vec);
}

//' Get the field named \code{data_has_been_set_} for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @param a boolean flag indicating if data has been set
// [[Rcpp::export(.LinOp__has_numerical_data)]]
bool LinOp__has_numerical_data(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  return ptr->has_numerical_data;
}

//' Get the \code{linOp_data_} field for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @return a LinOp XPtr
// [[Rcpp::export(.LinOp__get_linOp_data)]]
Rcpp::XPtr<LinOp> LinOp__get_linOp_data(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> linop_ptr(xp);
  Rcpp::XPtr<LinOp> ptr( linop_ptr->get_linOp_data());
  return ptr;
}

//' Set the \code{linOp_data_} field for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @param tree the LinOp data
// [[Rcpp::export(.LinOp__set_linOp_data)]]
void LinOp__set_linOp_data(SEXP xp, const Rcpp::XPtr<LinOp> tree) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> linop_ptr(xp);
  Rcpp::XPtr<LinOp> tree_ptr(tree);
  linop_ptr->set_linOp_data(tree_ptr);
}

//' Get the `data_ndim_` field for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @return the value of `data_ndim_`
// [[Rcpp::export(.LinOp__get_data_ndim)]]
int LinOp__get_data_ndim(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  return ptr->get_data_ndim();
}

//' Set the `data_ndim_` field of the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @param ndim the `data_ndim_` value 
// [[Rcpp::export(.LinOp__set_data_ndim)]]
void LinOp__set_data_ndim(SEXP xp, int) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  // Set the result
  ptr->set_data_ndim(ndim);
}

//' Is the data for LinOp sparse?
//'
//' @param xp the LinOp Object XPtr
//' @return TRUE or FALSE
// [[Rcpp::export(.LinOp__is_sparse)]]
bool LinOp__is_sparse(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  return ptr->is_sparse();
}

//' Get the `sparse_data_` field from the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @return a \link[Matrix]{dgCMatrix-class} object
// [[Rcpp::export(.LinOp__get_sparse_data)]]
Matrix LinOp__get_sparse_data(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  return ptr->get_sparse_data();
}

//' Get the field `dense_data_` for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @return a Matrix object
// [[Rcpp::export(.LinOp__get_dense_data)]]
Eigen::MatrixXd LinOp__get_dense_data(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  return ptr->get_dense_data();
}

//' Set the field \code{dense_data_} of the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @param denseMat a standard matrix object in R
// [[Rcpp::export(.LinOp__set_dense_data)]]
void LinOp__set_dense_data(SEXP xp, Eigen::Map<Eigen::MatrixXd> denseMat) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  // Set the result
  ptr->set_dense_data(denseMat);
}

//' Set the field `sparse_data_` of the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @param sparseMat a \link[Matrix]{dgCMatrix-class} object
// [[Rcpp::export(.LinOp__set_sparse_data)]]
void LinOp__set_sparse_data(SEXP xp, MappedSparseMatrix sparseMat) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  ptr->set_sparse_data(sparseMat);
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






