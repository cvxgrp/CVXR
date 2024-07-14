## A Constant vector of Linop Types
LINOP_TYPES <-
  c( ## from LinOp.hpp
    "VARIABLE",
    "PARAM",
    "PROMOTE",
    "MUL",
    "RMUL",
    "MUL_ELEM",
    "DIV",
    "SUM",
    "NEG",
    "INDEX",
    "TRANSPOSE",
    "SUM_ENTRIES",
    "TRACE",
    "RESHAPE",
    "DIAG_VEC",
    "DIAG_MAT",
    "UPPER_TRI",
    "CONV",
    "HSTACK",
    "VSTACK",
    "SCALAR_CONST",
    "DENSE_CONST",
    "SPARSE_CONST",
    "NO_OP",
    "KRON", ## for backwards compatibility; equivalent to KRON_R (which is preferred)
    "KRON_R",
    "KRON_L"
  )
names(LINOP_TYPES) <- LINOP_TYPES
