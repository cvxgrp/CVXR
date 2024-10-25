## CVXPY SOURCE: cvxpy/settings.py

##########################
#                        #
# CVXR Package Constants #
#                        #
##########################
# Constants for operators.
PLUS <- "+"
MINUS <- "-"
MUL <- "*"

# Prefix for default named variables.
VAR_PREFIX <- "var"
# Prefix for default named parameters.
PARAM_PREFIX <- "param"

# Used to overload ==.
NP_EQUAL_STR <- "equal"

# Constraint types.
EQ_CONSTR <- "=="
INEQ_CONSTR <- "<="

# Map of constraint types
# TODO: These should be defined in a solver model.
EQ_MAP <- "1"
LEQ_MAP <- "2"
SOC_MAP <- "3"
SOC_EW_MAP <- "4"
PSD_MAP <- "5"
EXP_MAP <- "6"
BOOL_MAP <- "7"
INT_MAP <- "8"

# Keys in the dictionary of cone dimensions.
# Cone dims are now defined in matrix stuffing modules rather than the solver module.

## EQ_DIM <- "f"   # SCS 3.0 changes this to "z"
EQ_DIM <- "z"
LEQ_DIM <- "l"
SOC_DIM <- "q"
PSD_DIM <- "s"
EXP_DIM <- "ep"

# Keys for non-convex constraints.
BOOL_IDS <- "bool_ids"
BOOL_IDX <- "bool_idx"
INT_IDS <- "int_ids"
INT_IDX <- "int_idx"

# Parametrized problem.
PARAM_PROB <- "param_prob"

# Keys for problem data dict.
C_KEY <- "c"
OFFSET <- "offset"
P_KEY <- "P"
Q_KEY <- "q"
A_KEY <- "A"
B_KEY <- "b"
G_KEY <- "G"
H_KEY <- "h"
F_KEY <- "F"
DIMS <- "dims"
BOOL_IDX <- "bool_vars_idx"
INT_IDX <- "int_vars_idx"

# Keys for curvature and sign
CONSTANT <- "CONSTANT"
AFFINE <- "AFFINE"
CONVEX <- "CONVEX"
CONCAVE <- "CONCAVE"
QUASILINEAR <- "QUASILINEAR"
QUASICONVEX <- "QUASICONVEX"
QUASICONCAVE <- "QUASICONCAVE"
LOG_LOG_CONSTANT <- "LOG-LOG CONSTANT"
LOG_LOG_AFFINE <- "LOG-LOG AFFINE"
LOG_LOG_CONVEX <- "LOG-LOG CONVEX"
LOG_LOG_CONCAVE <- "LOG-LOG CONCAVE"
ZERO <- "ZERO"
NONNEG <- "NONNEGATIVE"
NONPOS <- "NONPOSITIVE"
UNKNOWN <- "UNKNOWN"

# Canonicalization backends.
RUST_CANON_BACKEND <- "RUST"
CPP_CANON_BACKEND <- "CPP"
DEFAULT_CANON_BACKEND <- CPP_CANON_BACKEND

# Numerical tolerances.
EIGVAL_TOL <- 1e-10
PSD_NSD_PROJECTION_TOL <- 1e-8
GENERAL_PROJECTION_TOL <- 1e-10
SPARSE_PROJECTION_TOL <- 1e-10
ATOM_EVAL_TOL <- 1e-4

# Max number of nodes a reasonable expression should have (used for debugging).
MAX_NODES = 10000

# DPP is slow when total size of parameters exceed this threshold.
PARAM_THRESHOLD <- 1e4   # TODO: Should we reduce this?

# TODO: Can we set the number of threads to use during compilation like in CVXPY?
NUM_THREADS <- -1

