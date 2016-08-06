#include "LinOp.hpp"

//' Create a new LinOp object.
//'
//' @return an external ptr (Rcpp::XPtr) to a LinOp object instance.
// [[Rcpp::export]]
SEXP LinOp__new() {
  // create a pointer to an Uniform object and wrap it
  // as an external pointer
  Rcpp::XPtr<LinOp> ptr( new LinOp( ), true );

  // return the external pointer to the R side
  return ptr;
}

//' Get the field named \code{sparse_data} from the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @return a \link[Matrix]{dgCMatrix-class} object
// [[Rcpp::export]]
Eigen::SparseMatrix<double> LinOp__get_sparse_data(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  return ptr->sparse_data;
}

//' Set the field named \code{sparse_data} of the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @param sparseMat a \link[Matrix]{dgCMatrix-class} object
// [[Rcpp::export]]
void LinOp__set_sparse_data(SEXP xp, SEXP sparseMat) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  // Map to Eigen::SparseMatrix<double> which is Matrix courtesy of Utils.hpp typedef
  ptr->sparse_data = Rcpp::as<Eigen::SparseMatrix<double> >(sparseMat);

}

//' Get the field \code{dense_data} for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @return a MatrixXd object
// [[Rcpp::export]]
Eigen::MatrixXd LinOp__get_dense_data(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  return ptr->dense_data;
}

//' Set the field \code{dense_data} of the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @param denseMat a standard matrix object in R
// [[Rcpp::export]]
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
// [[Rcpp::export]]
std::vector<int>  LinOp__get_size(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  return ptr->size;
}

//' Set the field \code{size} of the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @param an integer vector object in R
// [[Rcpp::export]]
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
//' @param tree the LinOp Object XPtr to push
// [[Rcpp::export]]
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
// [[Rcpp::export]]
void LinOp__size_push_back(SEXP xp, int intVal) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  (ptr->size).push_back(intVal);
}

//' Set the field named \code{type} for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @param type an integer value
// [[Rcpp::export]]
void LinOp__set_type(SEXP xp, int typeValue) {
  OperatorType oType;
  int err = 0; // to signal error
  switch (typeValue) {
  case 0:
    oType = VARIABLE;
    break;
  case 1:
    oType = PROMOTE;
    break;
  case 2:
    oType = MUL;
    break;
  case 3:
    oType = RMUL;
    break;
  case 4:
    oType = MUL_ELEM;
    break;
  case 5:
    oType = DIV;
    break;
  case 6:
    oType = SUM;
    break;
  case 7:
    oType = NEG;
    break;
  case 8:
    oType = INDEX;
    break;
  case 9:
    oType = TRANSPOSE;
    break;
  case 10:
    oType = SUM_ENTRIES;
    break;
  case 11:
    oType = TRACE;
    break;
  case 12:
    oType = RESHAPE;
    break;
  case 13:
    oType = DIAG_VEC;
    break;
  case 14:
    oType = DIAG_MAT;
    break;
  case 15:
    oType = UPPER_TRI;
    break;
  case 16:
    oType = CONV;
    break;
  case 17:
    oType = HSTACK;
    break;
  case 18:
    oType = VSTACK;
    break;
  case 19:
    oType = SCALAR_CONST;
    break;
  case 20:
    oType = DENSE_CONST;
    break;
  case 21:
    oType = SPARSE_CONST;
    break;
  case 22:
    oType = NO_OP;
    break;
  case 23:
    oType = KRON;
    break;
  case 24:
    oType = EQ;     // equality constraint
    break;
  case 25:
    oType = LEQ;    // non-negative orthant
    break;
  case 26:
    oType = SOC;    // second-order cone
    break;
  case 27:
    oType = EXP;    // exponential cone
    break;
  case 28:
    oType = SDP;    // semi-definite cone
    break;
  default:
    err = 1;
    // std::cerr << "Error: linOp type invalid." << lin.type << std::endl;
    Rcpp::stop("Error: LinOp type invalid");
  }
  if (err < 1) {
    // grab the object as a XPtr (smart pointer)
    Rcpp::XPtr<LinOp> ptr(xp);
    ptr->type = oType;

#ifdef CVXCANON_DEBUG
    Rcpp::Rcout << "LinOp Id " << ptr->id << " type now " << ptr->type << std::endl;
#endif

  }
}

//' Get the field named \code{type} for the LinOp object
//'
//' @param xp the LinOp Object XPtr
//' @return an integer value for type
// [[Rcpp::export]]
int LinOp__get_type(SEXP xp) {

  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  int oType;

  switch (ptr->type) {
  case VARIABLE:
    oType = 0;
    break;
  case PROMOTE:
    oType = 1;
    break;
  case MUL:
    oType = 2;
    break;
  case RMUL:
    oType = 3;
    break;
  case MUL_ELEM:
    oType = 4;
    break;
  case DIV:
    oType = 5;
    break;
  case SUM:
    oType = 6;
    break;
  case NEG:
    oType = 7;
    break;
  case INDEX:
    oType = 8;
    break;
  case TRANSPOSE:
    oType = 9;
    break;
  case SUM_ENTRIES:
    oType = 10;
    break;
  case TRACE:
    oType = 11;
    break;
  case RESHAPE:
    oType = 12;
    break;
  case DIAG_VEC:
    oType = 13;
    break;
  case DIAG_MAT:
    oType = 14;
    break;
  case UPPER_TRI:
    oType = 15;
    break;
  case CONV:
    oType = 16;
    break;
  case HSTACK:
    oType = 17;
    break;
  case VSTACK:
    oType = 18;
    break;
  case SCALAR_CONST:
    oType = 19;
    break;
  case DENSE_CONST:
    oType = 20;
    break;
  case SPARSE_CONST:
    oType = 21;
    break;
  case NO_OP:
    oType = 22;
    break;
  case KRON:
    oType = 23;
    break;
  case EQ:
    oType = 24;
    break;
  case LEQ:
    oType = 25;
    break;
  case SOC:
    oType = 26;
    break;
  case EXP:
    oType = 27;
    break;
  case SDP:
    oType = 28;
    break;
  default:
    oType = -1;
    Rcpp::stop("Error: LinOp type invalid");
  }
  return oType;
}

//' Perform a push back operation on the \code{slice} field of LinOp
//'
//' @param xp the LinOp Object XPtr
//' @param intVec an intVector to push back
// [[Rcpp::export]]
void LinOp__slice_push_back(SEXP xp, std::vector<int> intVec) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);
  (ptr->slice).push_back(intVec);
}

//' Get the slice field of the LinOp Object
//'
//' @param xp the LinOp Object XPtr
//' @return the value of the slice field of the LinOp Object
// [[Rcpp::export]]
std::vector<std::vector<int> >  LinOp__get_slice(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  // Get the result
  return ptr->slice;
}

//' Get the slice field of the LinOp Object
//'
//' @param xp the LinOp Object XPtr
//' @return the value of the slice field of the LinOp Object
// [[Rcpp::export]]
std::string LinOp__get_id(SEXP xp) {
  // grab the object as a XPtr (smart pointer)
  Rcpp::XPtr<LinOp> ptr(xp);

  // Get the result
  return boost::lexical_cast<std::string>(ptr->id);
}






