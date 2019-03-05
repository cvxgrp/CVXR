#'
#' Sign Properties
#'
#' Determine if an expression is positive, negative, or zero.
#'
#' @param object An \linkS4class{Expression} object.
#' @return A logical value.
#' @examples
#' pos <- Constant(1)
#' neg <- Constant(-1)
#' zero <- Constant(0)
#' unknown <- Variable()
#'
#' is_zero(pos)
#' is_zero(-zero)
#' is_zero(unknown)
#' is_zero(pos + neg)
#'
#' is_nonneg(pos + zero)
#' is_nonneg(pos * neg)
#' is_nonneg(pos - neg)
#' is_nonneg(unknown)
#'
#' is_nonpos(-pos)
#' is_nonpos(pos + neg)
#' is_nonpos(neg * zero)
#' is_nonpos(neg - pos)
#' @name sign-methods
NULL

#' @rdname sign-methods
#' @export
setGeneric("is_zero", function(object) { standardGeneric("is_zero") })

#' @rdname sign-methods
#' @export
setGeneric("is_nonneg", function(object) { standardGeneric("is_nonneg") })

#' @rdname sign-methods
#' @export
setGeneric("is_nonpos", function(object) { standardGeneric("is_nonpos") })

#'
#' Complex Properties
#' 
#' Determine if an expression is real, imaginary, or complex.
#' 
#' @param object An \linkS4class{Expression} object.
#' @return A logical value.
#' @name complex-methods
NULL

#' @rdname complex-methods
#' @export
setGeneric("is_real", function(object) { standardGeneric("is_real") })

#' @rdname complex-methods
#' @export
setGeneric("is_imag", function(object) { standardGeneric("is_imag") })

#' @rdname complex-methods
#' @export
setGeneric("is_complex", function(object) { standardGeneric("is_complex") })

#'
#' Curvature of Expression
#'
#' The curvature of an expression.
#'
#' @param object An \linkS4class{Expression} object.
#' @return A string indicating the curvature of the expression, either "CONSTANT", "AFFINE", "CONVEX, "CONCAVE", or "UNKNOWN".
#' @examples
#' x <- Variable()
#' c <- Constant(5)
#'
#' curvature(c)
#' curvature(x)
#' curvature(x^2)
#' curvature(sqrt(x))
#' curvature(log(x^3) + sqrt(x))
#' @docType methods
#' @rdname curvature
#' @export
setGeneric("curvature", function(object) { standardGeneric("curvature") })

#'
#' Curvature Properties
#'
#' Determine if an expression is constant, affine, convex, concave, quadratic, piecewise linear (pwl), or quadratic/piecewise affine (qpwa).
#'
#' @param object An \linkS4class{Expression} object.
#' @return A logical value.
#' @examples
#' x <- Variable()
#' c <- Constant(5)
#'
#' is_constant(c)
#' is_constant(x)
#'
#' is_affine(c)
#' is_affine(x)
#' is_affine(x^2)
#'
#' is_convex(c)
#' is_convex(x)
#' is_convex(x^2)
#' is_convex(sqrt(x))
#'
#' is_concave(c)
#' is_concave(x)
#' is_concave(x^2)
#' is_concave(sqrt(x))
#'
#' is_quadratic(x^2)
#' is_quadratic(sqrt(x))
#'
#' is_pwl(c)
#' is_pwl(x)
#' is_pwl(x^2)
#' @name curvature-methods
NULL

#' @rdname curvature-methods
#' @export
setGeneric("is_constant", function(object) { standardGeneric("is_constant") })

#' @rdname curvature-methods
#' @export
setGeneric("is_affine", function(object) { standardGeneric("is_affine") })

#' @rdname curvature-methods
#' @export
setGeneric("is_convex", function(object) { standardGeneric("is_convex") })

#' @rdname curvature-methods
#' @export
setGeneric("is_concave", function(object) { standardGeneric("is_concave") })

#' @rdname curvature-methods
#' @export
setGeneric("is_quadratic", function(object) { standardGeneric("is_quadratic") })

#' @rdname curvature-methods
#' @export
setGeneric("is_pwl", function(object) { standardGeneric("is_pwl") })

#' @rdname curvature-methods
#' @export
setGeneric("is_qpwa", function(object) { standardGeneric("is_qpwa") })

#'
#' Log-Log Curvature of Expression
#'
#' The log-log curvature of an expression.
#'
#' @param object An \linkS4class{Expression} object.
#' @return A string indicating the log-log curvature of the expression, either "LOG_LOG_CONSTANT", "LOG_LOG_AFFINE", "LOG_LOG_CONVEX, "LOG_LOG_CONCAVE", or "UNKNOWN".
#' @docType methods
#' @rdname log_log_curvature
setGeneric("log_log_curvature", function(object) { standardGeneric("log_log_curvature") })

#'
#' Log-Log Curvature Properties
#'
#' Determine if an expression is log-log constant, log-log affine, log-log convex, or log-log concave.
#'
#' @param object An \linkS4class{Expression} object.
#' @return A logical value.
#' @name log_log_curvature-methods
NULL

#' @rdname log_log_curvature-methods
#' @export
setGeneric("is_log_log_constant", function(object) { standardGeneric("is_log_log_constant") })

#' @rdname log_log_curvature-methods
#' @export
setGeneric("is_log_log_affine", function(object) { standardGeneric("is_log_log_affine") })

#' @rdname log_log_curvature-methods
#' @export
setGeneric("is_log_log_convex", function(object) { standardGeneric("is_log_log_convex") })

#' @rdname log_log_curvature-methods
#' @export
setGeneric("is_log_log_concave", function(object) { standardGeneric("is_log_log_concave") })

#'
#' DCP Compliance
#'
#' Determine if a problem or expression complies with the disciplined convex programming rules.
#'
#' @param object A \linkS4class{Problem} or \linkS4class{Expression} object.
#' @return A logical value indicating whether the problem or expression is DCP compliant, i.e. no unknown curvatures.
#' @examples
#' x <- Variable()
#' prob <- Problem(Minimize(x^2), list(x >= 5))
#' is_dcp(prob)
#' solve(prob)
#' @docType methods
#' @rdname is_dcp
#' @export
setGeneric("is_dcp", function(object) { standardGeneric("is_dcp") })

#'
#' DGP Compliance
#' 
#' Determine if a problem or expression complies with the disciplined geometric programming rules.
#' 
#' @param object A \linkS4class{Problem} or \linkS4class{Expression} object.
#' @return A logical value indicating whether the problem or expression is DCP compliant, i.e. no unknown curvatures.
#' @examples 
#' x <- Variable()
#' y <- Variable()
#' prob <- Problem(Minimize(x*y), list(x <= 5, y >= -5))
#' is_dgp(prob)
#' solve(prob)
#' @docType methods
#' @rdname is_dgp
#' @export
setGeneric("is_dgp", function(object) { standardGeneric("is_dgp") })

#'
#' Size of Expression
#'
#' The size of an expression.
#'
#' @param object An \linkS4class{Expression} object.
#' @return A vector with two elements \code{c(row, col)} representing the dimensions of the expression.
#' @examples
#' x <- Variable()
#' y <- Variable(3)
#' z <- Variable(3,2)
#'
#' size(x)
#' size(y)
#' size(z)
#' size(x + y)
#' size(z - x)
#' @docType methods
#' @rdname size
#' @export
setGeneric("size", function(object) { standardGeneric("size") })

#'
#' Size Properties
#'
#' Determine if an expression is a scalar, vector, or matrix.
#'
#' @param object An \linkS4class{Expression} object.
#' @return A logical value.
#' @examples
#' x <- Variable()
#' y <- Variable(3)
#' z <- Variable(3,2)
#'
#' is_scalar(x)
#' is_scalar(y)
#' is_scalar(x + y)
#'
#' is_vector(x)
#' is_vector(y)
#' is_vector(2*z)
#'
#' is_matrix(x)
#' is_matrix(y)
#' is_matrix(z)
#' is_matrix(z - x)
#' @name size-methods
NULL

#' @rdname size-methods
#' @export
setGeneric("is_scalar", function(object) { standardGeneric("is_scalar") })

#' @rdname size-methods
#' @export
setGeneric("is_vector", function(object) { standardGeneric("is_vector") })

#' @rdname size-methods
#' @export
setGeneric("is_matrix", function(object) { standardGeneric("is_matrix") })

#'
#' Matrix Properties
#'
#' Determine if an expression is positive semidefinite, negative semidefinite, hermitian, and/or symmetric.
#' 
#' @param object An \linkS4class{Expression} object.
#' @return A logical value.
#' @name matrix_prop-methods
NULL

#' @rdname matrix_prop-methods
#' @export
setGeneric("is_psd", function(object) { standardGeneric("is_psd") })

#' @rdname matrix_prop-methods
#' @export
setGeneric("is_nsd", function(object) { standardGeneric("is_nsd") })

#' @rdname matrix_prop-methods
#' @export
setGeneric("is_hermitian", function(object) { standardGeneric("is_hermitian") })

#' @rdname matrix_prop-methods
#' @export
setGeneric("is_symmetric", function(object) { standardGeneric("is_symmetric") })

# The value of the objective given the solver primal value.
setGeneric("primal_to_result", function(object, result) { standardGeneric("primal_to_result") })

#'
#' Get or Set Value
#'
#' Get or set the value of a variable, parameter, expression, or problem.
#'
#' @param object A \linkS4class{Variable}, \linkS4class{Parameter}, \linkS4class{Expression}, or \linkS4class{Problem} object.
#' @param value A numeric scalar, vector, or matrix to assign to the object.
#' @return The numeric value of the variable, parameter, or expression. If any part of the mathematical object is unknown, return \code{NA}.
#' @examples
#' lambda <- Parameter()
#' value(lambda)
#'
#' value(lambda) <- 5
#' value(lambda)
#' @name value-methods
NULL

#' @rdname value-methods
#' @export
setGeneric("value", function(object) { standardGeneric("value") })

#' @rdname value-methods
#' @export
setGeneric("value<-", function(object, value) { standardGeneric("value<-") })

# Internal method for saving the value of an expression
setGeneric("save_value", function(object, value) { standardGeneric("save_value") })

#'
#' Get Expression Data
#'
#' Get information needed to reconstruct the expression aside from its arguments.
#'
#' @param object A \linkS4class{Expression} object.
#' @return A list containing data.
#' @docType methods
#' @rdname get_data
#' @export
setGeneric("get_data", function(object) { standardGeneric("get_data") })

#'
#' Variable, Parameter, or Expression Name
#'
#' The string representation of a variable, parameter, or expression.
#'
#' @param object A \linkS4class{Variable}, \linkS4class{Parameter}, or \linkS4class{Expression} object.
#' @return For \linkS4class{Variable} or \linkS4class{Parameter} objects, the value in the name slot. For \linkS4class{Expression} objects, a string indicating the nested atoms and their respective arguments.
#' @docType methods
#' @rdname name
#' @examples
#' x <- Variable()
#' y <- Variable(3, name = "yVar")
#'
#' name(x)
#' name(y)
#' @export
setGeneric("name", function(x) { standardGeneric("name") })

#'
#' Parts of an Expression Leaf
#'
#' List the variables, parameters, constants, or atoms in a canonical expression.
#'
#' @param object A \linkS4class{Leaf} object.
#' @return A list of \linkS4class{Variable}, \linkS4class{Parameter}, \linkS4class{Constant}, or \linkS4class{Atom} objects.
#' @examples
#' m <- 50
#' n <- 10
#' beta <- Variable(n)
#' y <- matrix(rnorm(m), nrow = m)
#' X <- matrix(rnorm(m*n), nrow = m, ncol = n)
#' lambda <- Parameter()
#'
#' expr <- sum_squares(y - X %*% beta) + lambda*p_norm(beta, 1)
#' variables(expr)
#' parameters(expr)
#' constants(expr)
#' lapply(constants(expr), function(c) { value(c) })
#' @name expression-parts
NULL

#' @rdname expression-parts
#' @export
setGeneric("variables", function(object) { standardGeneric("variables") })

#' @rdname expression-parts
#' @export
setGeneric("parameters", function(object) { standardGeneric("parameters") })

#' @rdname expression-parts
#' @export
setGeneric("constants", function(object) { standardGeneric("constants") })

#' @rdname expression-parts
#' @export
setGeneric("atoms", function(object) { standardGeneric("atoms") })

#'
#' Attributes of an Expression Leaf
#' 
#' Determine if an expression is positive or negative.
#' 
#' @param object A \linkS4class{Leaf} object.
#' @return A logical value.
#' @name leaf-attr
NULL

#' @rdname leaf-attr
#' @export
setGeneric("is_pos", function(object) { standardGeneric("is_pos") })

#' @rdname leaf-attr
#' @export
setGeneric("is_neg", function(object) { standardGeneric("is_neg") })

#'
#' Sub/Super-Gradient
#'
#' The (sub/super)-gradient of the expression with respect to each variable.
#' Matrix expressions are vectorized, so the gradient is a matrix. \code{NA} indicates variable values are unknown or outside the domain.
#'
#' @param object An \linkS4class{Expression} object.
#' @return A list mapping each variable to a sparse matrix.
#' @examples
#' x <- Variable(2, name = "x")
#' A <- Variable(2, 2, name = "A")
#'
#' value(x) <- c(-3,4)
#' expr <- p_norm(x, 2)
#' grad(expr)
#'
#' value(A) <- rbind(c(3,-4), c(4,3))
#' expr <- p_norm(A, 0.5)
#' grad(expr)
#'
#' value(A) <- cbind(c(1,2), c(-1,0))
#' expr <- abs(A)
#' grad(expr)
#' @docType methods
#' @rdname grad
#' @export
setGeneric("grad", function(object) { standardGeneric("grad") })

#'
#' Domain
#'
#' A list of constraints describing the closure of the region where the expression is finite.
#'
#' @param object An \linkS4class{Expression} object.
#' @return A list of \linkS4class{Constraint} objects.
#' @examples
#' a <- Variable(name = "a")
#' dom <- domain(p_norm(a, -0.5))
#' prob <- Problem(Minimize(a), dom)
#' result <- solve(prob)
#' result$value
#'
#' b <- Variable()
#' dom <- domain(kl_div(a, b))
#' result <- solve(Problem(Minimize(a + b), dom))
#' result$getValue(a)
#' result$getValue(b)
#'
#' A <- Variable(2, 2, name = "A")
#' dom <- domain(lambda_max(A))
#' A0 <- rbind(c(1,2), c(3,4))
#' result <- solve(Problem(Minimize(norm2(A - A0)), dom))
#' result$getValue(A)
#'
#' dom <- domain(log_det(A + diag(rep(1,2))))
#' prob <- Problem(Minimize(sum(diag(A))), dom)
#' result <- solve(prob, solver = "SCS")
#' result$value
#' @docType methods
#' @rdname domain
#' @export
setGeneric("domain", function(object) { standardGeneric("domain") })

#'
#' Project Value
#' 
#' Project a value onto the attribute set of a \linkS4class{Leaf}.
#' A sensible idiom is \code{value(leaf) = project(leaf, val)}.
#' 
#' @param object A \linkS4class{Leaf} object.
#' @param val The assigned value.
#' @return The value rounded to the attribute type.
#' @docType methods
#' @rdname project-methods
NULL

#' @rdname project-methods
#' @export
setGeneric("project", function(object, value) { standardGeneric("project") })

#' @rdname project-methods
#' @export
setGeneric("project_and_assign", function(object, value) { standardGeneric("project_and_assign") })

#'
#' Validate Value
#'
#' Check that the value satisfies a \linkS4class{Leaf}'s symbolic attributes.
#'
#' @param object A \linkS4class{Leaf} object.
#' @param val The assigned value.
#' @return The value converted to proper matrix type.
#' @docType methods
#' @rdname validate_val
setGeneric("validate_val", function(object, val) { standardGeneric("validate_val") })

#'
#' Canonicalize
#'
#' Computes the graph implementation of a canonical expression.
#'
#' @param object A \linkS4class{Canonical} object.
#' @return A list of \code{list(affine expression, list(constraints))}.
#' @docType methods
#' @name canonicalize
NULL

#' @rdname canonicalize
#' @export
setGeneric("canonicalize", function(object) { standardGeneric("canonicalize") })

#' @rdname canonicalize
setGeneric("canonical_form", function(object) { standardGeneric("canonical_form") })

#
# Gradient of an Atom
#
# The (sub/super) gradient of the atom with respect to each argument. Matrix expressions are vectorized, so the gradient is a matrix.
#
# @param object An \linkS4class{Atom} object.
# @param values A list of numeric values for the arguments.
# @return A list of sparse matrices or \code{NA}.
setGeneric(".grad", function(object, values) { standardGeneric(".grad") })

#
# Domain of an Atom
#
# The constraints describing the domain of the atom.
#
# @param object An \linkS4class{Atom} object.
# @return A list of \linkS4class{Constraint} objects.
setGeneric(".domain", function(object) { standardGeneric(".domain") })

#
# Gradient of an AxisAtom
#
# The (sub/super) gradient of the atom with respect to each argument. Matrix expressions are vectorized, so the gradient is a matrix. Takes the axis into account.
#
# @param values A list of numeric values for the arguments.
# @return A list of sparse matrices or \code{NA}.
setGeneric(".axis_grad", function(object, values) { standardGeneric(".axis_grad") })

#
# Column Gradient of an Atom
#
# The (sub/super) gradient of the atom with respect to a column argument. Matrix expressions are vectorized, so the gradient is a matrix.
# @param value A numeric value for a column.
# @return A sparse matrix or \code{NA}.
setGeneric(".column_grad", function(object, value) { standardGeneric(".column_grad") })

# Positive definite inequalities
#' @rdname PSDConstraint-class
#' @export
setGeneric("%>>%", function(e1, e2) { standardGeneric("%>>%") })

#' @rdname PSDConstraint-class
#' @export
setGeneric("%<<%", function(e1, e2) { standardGeneric("%<<%") })

#'
#' Atom Dimensions
#'
#' Determine the dimensions of an atom based on its arguments.
#'
#' @param object A \linkS4class{Atom} object.
#' @return A numeric vector \code{c(row, col)} indicating the dimensions of the atom.
#' @rdname dim_from_args
setGeneric("dim_from_args", function(object) { standardGeneric("dim_from_args") })

#'
#' Atom Sign
#'
#' Determine the sign of an atom based on its arguments.
#'
#' @param object An \linkS4class{Atom} object.
#' @return A logical vector \code{c(is positive, is negative)} indicating the sign of the atom.
#' @rdname sign_from_args
setGeneric("sign_from_args", function(object) { standardGeneric("sign_from_args") })

#'
#' Validate Arguments
#'
#' Validate an atom's arguments, returning an error if any are invalid.
#'
#' @param object An \linkS4class{Atom} object.
#' @docType methods
#' @rdname validate_args
setGeneric("validate_args", function(object) { standardGeneric("validate_args") })

#'
#' Numeric Value of Atom
#'
#' Returns the numeric value of the atom evaluated on the specified arguments.
#'
#' @param object An \linkS4class{Atom} object.
#' @param values A list of arguments to the atom.
#' @return A numeric scalar, vector, or matrix.
#' @docType methods
#' @rdname to_numeric
setGeneric("to_numeric", function(object, values) { standardGeneric("to_numeric") })

#'
#' Curvature of an Atom
#'
#' Determine if an atom is convex, concave, or affine.
#'
#' @param object A \linkS4class{Atom} object.
#' @return A logical value.
#' @examples
#' x <- Variable()
#'
#' is_atom_convex(x^2)
#' is_atom_convex(sqrt(x))
#' is_atom_convex(log(x))
#'
#' is_atom_concave(-abs(x))
#' is_atom_concave(x^2)
#' is_atom_concave(sqrt(x))
#'
#' is_atom_affine(2*x)
#' is_atom_affine(x^2)
#' @name curvature-atom
NULL

#' @rdname curvature-atom
#' @export
setGeneric("is_atom_convex", function(object) { standardGeneric("is_atom_convex") })

#' @rdname curvature-atom
#' @export
setGeneric("is_atom_concave", function(object) { standardGeneric("is_atom_concave") })

#' @rdname curvature-atom
#' @export
setGeneric("is_atom_affine", function(object) { standardGeneric("is_atom_affine") })

#'
#' Log-Log Curvature of an Atom
#'
#' Determine if an atom is log-log convex, concave, or affine.
#'
#' @param object A \linkS4class{Atom} object.
#' @return A logical value.
#' @name log_log_curvature-atom
NULL

#' @rdname log_log_curvature-atom
#' @export
setGeneric("is_atom_log_log_convex", function(object) { standardGeneric("is_atom_log_log_convex") })

#' @rdname log_log_curvature-atom
#' @export
setGeneric("is_atom_log_log_concave", function(object) { standardGeneric("is_atom_log_log_concave") })

#' @rdname log_log_curvature-atom
#' @export
setGeneric("is_atom_log_log_affine", function(object) { standardGeneric("is_atom_log_log_affine") })

#'
#' Curvature of Composition
#'
#' Determine whether a composition is non-decreasing or non-increasing in an index.
#'
#' @param object A \linkS4class{Atom} object.
#' @param idx An index into the atom.
#' @return A logical value.
#' @examples
#' x <- Variable()
#' is_incr(log(x), 1)
#' is_incr(x^2, 1)
#' is_decr(min(x), 1)
#' is_decr(abs(x), 1)
#' @name curvature-comp
NULL

#' @rdname curvature-comp
#' @export
setGeneric("is_incr", function(object, idx) { standardGeneric("is_incr") })

#' @rdname curvature-comp
#' @export
setGeneric("is_decr", function(object, idx) { standardGeneric("is_decr") })

#'
#' Graph Implementation
#'
#' Reduces the atom to an affine expression and list of constraints.
#'
#' @param object An \linkS4class{Expression} object.
#' @param arg_objs A list of linear expressions for each argument.
#' @param size A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @return A list of \code{list(LinOp for objective, list of constraints)}, where LinOp is a list representing the linear operator.
#' @docType methods
#' @rdname graph_implementation
setGeneric("graph_implementation", function(object, arg_objs, dim, data) { standardGeneric("graph_implementation") })

#'
#' Identification Number
#'
#' A unique identification number used internally to keep track of variables and constraints. Should not be modified by the user.
#'
#' @param object A \linkS4class{Variable} or \linkS4class{Constraint} object.
#' @return A non-negative integer identifier.
#' @seealso \code{\link[CVXR]{get_id}} \code{\link[CVXR]{setIdCounter}}
#' @examples
#' x <- Variable()
#' constr <- (x >= 5)
#' id(x)
#' id(constr)
#' @docType methods
#' @rdname id
#' @export
setGeneric("id", function(object) { standardGeneric("id") })

#'
#' Constraint Residual
#'
#' The residual expression of a constraint, i.e. the amount by which it is violated, and the value of that violation.
#' For instance, if our constraint is \eqn{g(x) \leq 0}, the residual is \eqn{max(g(x), 0)} applied elementwise.
#'
#' @param object A \linkS4class{Constraint} object.
#' @return A \linkS4class{Expression} representing the residual, or the value of this expression.
#' @docType methods
#' @name residual-methods
NULL

#' @rdname residual-methods
setGeneric("residual", function(object) { standardGeneric("residual") })

#' @rdname residual-methods
setGeneric("violation", function(object) { standardGeneric("violation") })

#'
#' Second-Order Cone Methods
#'
#' The number of elementwise cones or a list of the sizes of the elementwise cones.
#'
#' @param object An \linkS4class{SOCAxis} object.
#' @return The number of cones, or the size of a cone.
#' @docType methods
#' @name cone-methods
NULL

#' @rdname cone-methods
setGeneric("num_cones", function(object) { standardGeneric("num_cones") })

#' @rdname cone-methods
setGeneric("cone_sizes", function(object) { standardGeneric("cone_sizes") })

#
# Get and Set Dual Value
#
# Get and set the value of the dual variable in a constraint.
#
# @param object A \linkS4class{Constraint} object.
# @param value A numeric scalar, vector, or matrix to assign to the object.
# @return The numeric value of the dual variable. Defaults to \code{NA} if unknown.
# @rdname dual_value-methods
# NULL

# @rdname dual_value-methods
# @export
setGeneric("dual_value", function(object) { standardGeneric("dual_value") })

# @rdname dual_value-methods
# @export
setGeneric("dual_value<-", function(object, value) { standardGeneric("dual_value<-") })

#'
#' Format Constraints
#'
#' Format constraints for the solver.
#'
#' @param object A \linkS4class{Constraint} object.
#' @param eq_constr A list of the equality constraints in the canonical problem.
#' @param leq_constr A list of the inequality constraints in the canonical problem.
#' @param dims A list with the dimensions of the conic constraints.
#' @param solver A string representing the solver to be called.
#' @return A list containing equality constraints, inequality constraints, and dimensions.
#' @rdname format_constr
setGeneric("format_constr", function(object, eq_constr, leq_constr, dims, solver) { standardGeneric("format_constr") })

# Constraint generic methods
setGeneric("constr_type", function(object) { standardGeneric("constr_type") })
setGeneric("constr_id", function(object) { standardGeneric("constr_id") })

# Nonlinear constraint generic methods
setGeneric("block_add", function(object, mat, block, vert_offset, horiz_offset, rows, cols, vert_step, horiz_step) { standardGeneric("block_add") })
setGeneric("place_x0", function(object, big_x, var_offsets) { standardGeneric("place_x0") })
setGeneric("place_Df", function(object, big_Df, Df, var_offsets, vert_offset) { standardGeneric("place_Df") })
setGeneric("place_H", function(object, big_H, H, var_offsets) { standardGeneric("place_H") })
setGeneric("extract_variables", function(object, x, var_offsets) { standardGeneric("extract_variables") })

# Problem generic methods
#'
#' Parts of a Problem
#'
#' Get and set the objective, constraints, or size metrics (get only) of a problem.
#'
#' @param object A \linkS4class{Problem} object.
#' @param value The value to assign to the slot.
#' @return For getter functions, the requested slot of the object.
#' x <- Variable()
#' prob <- Problem(Minimize(x^2), list(x >= 5))
#' objective(prob)
#' constraints(prob)
#' size_metrics(prob)
#'
#' objective(prob) <- Maximize(sqrt(x))
#' constraints(prob) <- list(x <= 10)
#' objective(prob)
#' constraints(prob)
#' @name problem-parts
NULL

#' @rdname problem-parts
#' @export
setGeneric("objective", function(object) { standardGeneric("objective") })

#' @rdname problem-parts
#' @export
setGeneric("objective<-", function(object, value) { standardGeneric("objective<-") })

#' @rdname problem-parts
#' @export
setGeneric("constraints", function(object) { standardGeneric("constraints") })

#' @rdname problem-parts
#' @export
setGeneric("constraints<-", function(object, value) { standardGeneric("constraints<-") })

#' @rdname problem-parts
#' @export
setGeneric("size_metrics", function(object) { standardGeneric("size_metrics") })

setGeneric("status", function(object) { standardGeneric("status") })
setGeneric("status<-", function(object, value) { standardGeneric("status<-") })

setGeneric("solver_stats", function(object) { standardGeneric("solver_stats") })
setGeneric("solver_stats<-", function(object, value) { standardGeneric("solver_stats<-") })

#'
#' Get Problem Data
#'
#' Get the problem data used in the call to the solver.
#'
#' @param object A \linkS4class{Problem} object.
#' @param solver A string indicating the solver that the problem data is for. Call \code{installed_solvers()} to see all available.
#' @return A list of arguments for the solver.
#' @examples
#' a <- Variable(name = "a")
#' data <- get_problem_data(Problem(Maximize(exp(a) + 2)), "SCS")
#' data[["dims"]]
#' data[["c"]]
#' data[["A"]]
#'
#' x <- Variable(2, name = "x")
#' data <- get_problem_data(Problem(Minimize(p_norm(x) + 3)), "ECOS")
#' data[["dims"]]
#' data[["c"]]
#' data[["A"]]
#' data[["G"]]
#' @rdname get_problem_data
#' @export
setGeneric("get_problem_data", function(object, solver, gp) { standardGeneric("get_problem_data") })

#'
#' Solve a DCP Problem
#'
#' Solve a DCP compliant optimization problem.
#'
#' @param object,a A \linkS4class{Problem} object.
#' @param solver,b (Optional) A string indicating the solver to use. Defaults to "ECOS".
#' @param ignore_dcp (Optional) A logical value indicating whether to override the DCP check for a problem.
#' @param warm_start (Optional) A logical value indicating whether the previous solver result should be used to warm start.
#' @param verbose (Optional) A logical value indicating whether to print additional solver output.
#' @param parallel (Optional) A logical value indicating whether to solve in parallel if the problem is separable.
#' @param ... Additional options that will be passed to the specific solver. In general, these options will override any default settings imposed by CVXR.
#' @return A list containing the solution to the problem:
#' \describe{
#'    \item{\code{status}}{The status of the solution. Can be "optimal", "optimal_inaccurate", "infeasible", "infeasible_inaccurate", "unbounded", "unbounded_inaccurate", or "solver_error".}
#'    \item{\code{value}}{The optimal value of the objective function.}
#'    \item{\code{solver}}{The name of the solver.}
#'    \item{\code{solve_time}}{The time (in seconds) it took for the solver to solve the problem.}
#'    \item{\code{setup_time}}{The time (in seconds) it took for the solver to set up the problem.}
#'    \item{\code{num_iters}}{The number of iterations the solver had to go through to find a solution.}
#'    \item{\code{getValue}}{A function that takes a \linkS4class{Variable} object and retrieves its primal value.}
#'    \item{\code{getDualValue}}{A function that takes a \linkS4class{Constraint} object and retrieves its dual value(s).}
#' }
#' @examples
#' a <- Variable(name = "a")
#' prob <- Problem(Minimize(norm_inf(a)), list(a >= 2))
#' result <- psolve(prob, solver = "ECOS", verbose = TRUE)
#' result$status
#' result$value
#' result$getValue(a)
#' result$getDualValue(constraints(prob)[[1]])
#' @docType methods
#' @aliases psolve solve
#' @rdname psolve
#' @export
setGeneric("psolve", function(object, solver, ignore_dcp = FALSE, warm_start = FALSE, verbose = FALSE, parallel = FALSE, ...) { standardGeneric("psolve") })

#'
#' Is Problem a QP?
#'
#' Determine if a problem is a quadratic program.
#'
#' @param object A \linkS4class{Problem} object.
#' @return A logical value indicating whether the problem is a quadratic program.
#' @docType methods
#' @rdname is_qp
#' @export
setGeneric("is_qp", function(object) { standardGeneric("is_qp") })

#'
#' Is Problem Mixed Integer?
#' 
#' Determine if a problem is a mixed-integer program.
#' 
#' @param object A \linkS4class{Problem} object.
#' @return A logical value indicating whether the problem is a mixed-integer program
#' @docType methods
#' @rdname is_mixed_integer
#' @export
setGeneric("is_mixed_integer", function(object) { standardGeneric("is_mixed_integer") })

#'
#' Parse output from a solver and updates problem state
#'
#' Updates problem status, problem value, and primal and dual variable values
#'
#' @param object A \linkS4class{Problem} object.
#' @param solver A character string specifying the solver such as "ECOS", "SCS" etc.
#' @param results_dict the solver output
#' @return A list containing the solution to the problem:
#' \describe{
#'    \item{\code{status}}{The status of the solution. Can be "optimal", "optimal_inaccurate", "infeasible", "infeasible_inaccurate", "unbounded", "unbounded_inaccurate", or "solver_error".}
#'    \item{\code{value}}{The optimal value of the objective function.}
#'    \item{\code{solver}}{The name of the solver.}
#'    \item{\code{solve_time}}{The time (in seconds) it took for the solver to solve the problem.}
#'    \item{\code{setup_time}}{The time (in seconds) it took for the solver to set up the problem.}
#'    \item{\code{num_iters}}{The number of iterations the solver had to go through to find a solution.}
#'    \item{\code{getValue}}{A function that takes a \linkS4class{Variable} object and retrieves its primal value.}
#'    \item{\code{getDualValue}}{A function that takes a \linkS4class{Constraint} object and retrieves its dual value(s).}
#' }
#' @examples
#' \dontrun{
#' x <- Variable(2)
#' obj <- Minimize(x[1] + cvxr_norm(x, 1))
#' constraints <- list(x >= 2)
#' prob1 <- Problem(obj, constraints)
#' # Solve with ECOS.
#' ecos_data <- get_problem_data(prob1, "ECOS")
#' # Call ECOS solver interface directly
#' ecos_output <- ECOSolveR::ECOS_csolve(
#'                            c = ecos_data[["c"]],
#'                            G = ecos_data[["G"]],
#'                            h = ecos_data[["h"]],
#'                            dims = ecos_data[["dims"]],
#'                            A = ecos_data[["A"]],
#'                            b = ecos_data[["b"]]
#'                          )
#' # Unpack raw solver output.
#' res1 <- unpack_results(prob1, "ECOS", ecos_output)
#' # Without DCP validation (so be sure of your math), above is equivalent to:
#' # res1 <- solve(prob1, solver = "ECOS")
#' X <- Semidef(2)
#' Fmat <- rbind(c(1,0), c(0,-1))
#' obj <- Minimize(sum_squares(X - Fmat))
#' prob2 <- Problem(obj)
#' scs_data <- get_problem_data(prob2, "SCS")
#' scs_output <- scs::scs(
#'                       A = scs_data[['A']],
#'                       b = scs_data[['b']],
#'                       obj = scs_data[['c']],
#'                       cone = scs_data[['dims']]
#'                   )
#' res2 <- unpack_results(prob2, "SCS", scs_output)
#' # Without DCP validation (so be sure of your math), above is equivalent to:
#' # res2 <- solve(prob2, solver = "SCS")
#' }
#' @docType methods
#' @rdname unpack_results
#' @export
setGeneric("unpack_results", function(object, solver, results_dict) { standardGeneric("unpack_results") })
setGeneric(".handle_no_solution", function(object, status) { standardGeneric(".handle_no_solution") })
setGeneric("Problem.save_values", function(object, result_vec, objstore, offset_map) { standardGeneric("Problem.save_values") })
setGeneric(".save_dual_values", function(object, result_vec, constraints, constr_types) { standardGeneric(".save_dual_values") })
setGeneric(".update_problem_state", function(object, results_dict, sym_data, solver) { standardGeneric(".update_problem_state") })

# Problem data generic methods
setGeneric("get_objective", function(object) { standardGeneric("get_objective") })
setGeneric("get_eq_constr", function(object) { standardGeneric("get_eq_constr") })
setGeneric("get_ineq_constr", function(object) { standardGeneric("get_ineq_constr") })
setGeneric("get_nonlin_constr", function(object) { standardGeneric("get_nonlin_constr") })

# Solver generic methods
#'
#' Import Solver
#'
#' Import the R library that interfaces with the specified solver.
#'
#' @param solver A \linkS4class{Solver} object.
#' @examples
#' import_solver(ECOS())
#' import_solver(SCS())
#' @rdname import_solver
#' @docType methods
#' @export
setGeneric("import_solver", function(solver) { standardGeneric("import_solver") })
setGeneric("is_installed", function(solver) { standardGeneric("is_installed") })

#
# Choose MOSEK Solution
#
# Chooses between the basic and interior point solution from MOSEK. Solutions are ranked optimal > near_optimal > anything else > None.
# As long as interior solution is not worse, take it (for backward compatibility).
#
# @param solver A \linkS4class{MOSEK} object.
# @param results_dict A list of the results returned by the solver.
# @return A list containing the preferred solution (\code{solist}) and status of the preferred solution (\code{solsta}).
# @rdname choose_solution
# setGeneric("choose_solution", function(solver, results_dict) { standardGeneric("choose_solution") })

setGeneric("nonlin_constr", function(solver) { standardGeneric("nonlin_constr") })

#'
#' Validate Solver
#'
#' Raises an exception if the solver cannot solve the problem.
#'
#' @param solver A \linkS4class{Solver} object.
#' @param problem A \linkS4class{Problem} object.
#' @docType methods
#' @rdname validate_solver
setGeneric("validate_solver", function(solver, problem) { standardGeneric("validate_solver") })

#
# Validate Cache
#
# Clears the cache if the objective or constraints changed.
#
# @param solver A \linkS4class{Solver} object.
# @param objective A list representing the canonicalized objective.
# @param constraints A list of canonicalized constraints.
# @param cached_data A list mapping solver name to cached problem data.
# @return The updated \code{cached_data}.
setGeneric("validate_cache", function(solver, objective, constraints, cached_data) { standardGeneric("validate_cache") })

#
# Get Symbolic Data
#
# Returns the symbolic data for the problem.
#
# @param solver A \linkS4class{Solver} object.
# @param objective A list representing the canonicalized objective.
# @param constraints A list of canonicalized constraints.
# @param cached_data A list mapping solver name to cached problem data.
# @return A \linkS4class{SymData} object holding the symbolic data for the problem.
setGeneric("get_sym_data", function(solver, objective, constraints, cached_data) { standardGeneric("get_sym_data") })

#
# Get Matrix Data
#
# Returns the numeric data for the problem.
#
# @param solver A \linkS4class{Solver} object.
# @param objective A list representing the canonicalized objective.
# @param constraints A list of canonicalized constraints.
# @param cached_data A list mapping solver name to cached problem data.
# @return A \linkS4class{SymData} object holding the symbolic data for the problem.
setGeneric("get_matrix_data", function(solver, objective, constraints, cached_data) { standardGeneric("get_matrix_data") })

#
# Solver: Get Problem Data
#
# Returns the argument for the call to the solver.
#
# @param solver A \linkS4class{Solver} object.
# @param objective A list representing the canonicalized objective.
# @param constraints A list of canonicalized constraints.
# @param cached_data A list mapping solver name to cached problem data.
# @return A list of the arguments needed for the solver.
setGeneric("Solver.get_problem_data", function(solver, objective, constraints, cached_data) { standardGeneric("Solver.get_problem_data") })

#
# Extract Constraints
#
# Extracts the equality, inequality, and nonlinear constraints.
#
# @param solver A \linkS4class{Solver} object.
# @param constr_map A list of canonicalized constraints.
# @return A list of equality, inequality, and nonlinear constraints.
setGeneric("split_constr", function(solver, constr_map) { standardGeneric("split_constr") })

#'
#' Call to Solver
#'
#' Returns the result of the call to the solver.
#'
#' @param solver A \linkS4class{Solver} object.
#' @param objective A list representing the canonicalized objective.
#' @param constraints A list of canonicalized constraints.
#' @param cached_data A list mapping solver name to cached problem data.
#' @param warm_start A logical value indicating whether the previous solver result should be used to warm start.
#' @param verbose A logical value indicating whether to print solver output.
#' @param ... Additional arguments to the solver.
#' @return A list containing the status, optimal value, primal variable, and dual variables for the equality and inequality constraints.
#' @docType methods
#' @aliases Solver.solve
#' @rdname Solver-solve
setGeneric("Solver.solve", function(solver, objective, constraints, cached_data, warm_start, verbose, ...) { standardGeneric("Solver.solve") })

#'
#' Format Solver Results
#'
#' Converts the solver output into standard form.
#'
#' @param solver A \linkS4class{Solver} object.
#' @param results_dict A list containing the solver output.
#' @param data A list containing information about the problem.
#' @param cached_data A list mapping solver name to cached problem data.
#' @return A list containing the solver output in standard form.
#' @docType methods
#' @rdname format_results
setGeneric("format_results", function(solver, results_dict, data, cached_data) { standardGeneric("format_results") })

#'
#' Solver Capabilities
#'
#' Determine if a solver is capable of solving a linear program (LP), second-order cone program (SOCP), positive semidefinite program (PSD), exponential cone program (EXP), or mixed-integer program (MIP).
#'
#' @param solver A \linkS4class{Solver} object.
#' @return A logical value.
#' @examples
#' lp_capable(ECOS())
#' socp_capable(ECOS())
#' psd_capable(ECOS())
#' exp_capable(ECOS())
#' mip_capable(ECOS())
#' @name Solver-capable
NULL

#' @rdname Solver-capable
#' @export
setGeneric("lp_capable", function(solver) { standardGeneric("lp_capable") })

#' @rdname Solver-capable
#' @export
setGeneric("socp_capable", function(solver) { standardGeneric("socp_capable") })

#' @rdname Solver-capable
#' @export
setGeneric("psd_capable", function(solver) { standardGeneric("psd_capable") })

#' @rdname Solver-capable
#' @export
setGeneric("exp_capable", function(solver) { standardGeneric("exp_capable") })

#' @rdname Solver-capable
#' @export
setGeneric("mip_capable", function(solver) { standardGeneric("mip_capable") })

# Map of solver status code to CVXR status.
setGeneric("status_map", function(solver, status) { standardGeneric("status_map") })

# Map of CBC MIP/LP status to CVXR status.
setGeneric("status_map_mip", function(solver, status) { standardGeneric("status_map_mip") })
setGeneric("status_map_lp", function(solver, status) { standardGeneric("status_map_lp") })

#'
#' Reduction Acceptance
#' 
#' Determine whether the reduction accepts a problem.
#' 
#' @param object A \linkS4class{Reduction} object.
#' @param problem A \linkS4class{Problem} to check.
#' @return A logical value indicating whether the reduction can be applied.
#' @docType methods
#' @rdname accepts
setGeneric("accepts", function(object, problem) { standardGeneric("accepts") })

#'
#' Reduce a Problem
#'
#' Reduces the owned problem to an equivalent problem.
#' 
#' @param object A \linkS4class{Reduction} object.
#' @return An equivalent problem, encoded either as a \linkS4class{Problem} object or a list.
#' @docType methods
#' @rdname reduce
setGeneric("reduce", function(object) { standardGeneric("reduce") })

#'
#' Retrieve Solution
#' 
#' Retrieves a solution to the owned problem.
#' 
#' @param object A \linkS4class{Reduction} object.
#' @return A \linkS4class{Solution} to the problem emitted by \code{\link{reduce}}.
#' @docType methods
#' @rdname retrieve
setGeneric("retrieve", function(object, solution) { standardGeneric("retrieve") })

#'
#' Perform Reduction
#' 
#' Performs the reduction on a problem and returns an equivalent problem.
#' 
#' @param object A \linkS4class{Reduction} object.
#' @param problem A \linkS4class{Problem} on which the reduction will be performed.
#' @return A list containing
#' \itemize{
#'    \item{"problem"}{A \linkS4class{Problem} or list representing the equivalent problem.}
#'    \item{"inverse_data"}{A \linkS4class{InverseData} or list containing the data needed to invert this particular reduction.}
#' }
#' @docType methods
#' @rdname perform
setGeneric("perform", function(object, problem) { standardGeneric("perform") })

#'
#' Return Original Solution
#' 
#' Returns a solution to the original problem given the inverse data.
#' 
#' @param object A \linkS4class{Reduction} object.
#' @param solution A \linkS4class{Solution} to a problem that generated \code{inverse_data}.
#' @param inverse_data A \linkS4class{InverseData} object encoding the original problem.
#' @return A \linkS4class{Solution} to the original problem.
#' @docType methods
#' @rdname invert
setGeneric("invert", function(object, solution, inverse_data) { standardGeneric("invert") })
setGeneric("mosek_invert", function(object, results, inverse_data) { standardGeneric("mosek_invert") })

## Version 1.0 edits

## Start of newly added generics
setGeneric("expr", function(object) { standardGeneric("expr") })
setGeneric("ndim", function(object) { standardGeneric("ndim") })
setGeneric("flatten", function(object) { standardGeneric("flatten") })
setGeneric("var_id", function(object) { standardGeneric("var_id") })
setGeneric("get_attr_str", function(object) { standardGeneric("get_attr_str") })
setGeneric("get_var_offsets", function(object, variables) { standardGeneric("get_var_offsets") })
setGeneric("allow_complex", function(object) { standardGeneric("allow_complex") })
setGeneric("canonicalize_tree", function(object, expr) { standardGeneric("canonicalize_tree") })
setGeneric("canonicalize_expr", function(object, expr, args) { standardGeneric("canonicalize_expr") })
setGeneric("stuffed_objective", function(object, problem, inverse_data) { standardGeneric("stuffed_objective") })
setGeneric("cone_stuffed_objective", function(object, problem, extractor) { standardGeneric("cone_stuffed_objective") })
setGeneric("prepend", function(object, chain) { standardGeneric("prepend") })
setGeneric("group_coeff_offset", function(object, problem, constraints, exp_cone_order) { standardGeneric("group_coeff_offset") })
setGeneric("construct_intermediate_chain", function(problem, candidates, gp) { standardGeneric("construct_intermediate_chain") })
setGeneric("solve_via_data", function(object, data, warm_start, verbose, solver_opts, solver_cache) { standardGeneric("solve_via_data") })
setGeneric("unpack", function(object, solution) { standardGeneric("unpack") })
setGeneric("unpack_results", function(object, solution, chain, inverse_data) { standardGeneric("unpack_results") })
setGeneric("constr_value", function(object, tolerance) { standardGeneric("constr_value") })
setGeneric(".construct_dual_variables", function(object, args) { standardGeneric(".construct_dual_variables") })
setGeneric("reduction_solve", function(object, problem, warm_start, verbose, solver_opts) { standardGeneric("reduction_solve") })
setGeneric("reduction_solve_via_data", function(object, problem, data, warm_start, verbose, solver_opts) { standardGeneric("reduction_solve_via_data") })
setGeneric("reduction_format_constr", function(object, problem, constr, exp_cone_order) { standardGeneric("reduction_format_constr") })
setGeneric("block_format", function(object, problem, constraints, exp_cone_order) { standardGeneric("block_format") })
setGeneric("suitable", function(solver, problem) { standardGeneric("suitable") })
setGeneric("ls_solve", function(object, objective, constraints, cached_data, warm_start, verbose, solver_opts) { standardGeneric("ls_solve") })
setGeneric("supported_constraints", function(solver) { standardGeneric("supported_constraints") })
setGeneric("requires_constr", function(solver) { standardGeneric("requires_constr") })
setGeneric("get_coeffs", function(object, expr) { standardGeneric("get_coeffs") })
setGeneric("constant", function(object, expr) { standardGeneric("constant") })
setGeneric("affine", function(object, expr) { standardGeneric("affine") })
setGeneric("extract_quadratic_coeffs", function(object, affine_expr, quad_forms) { standardGeneric("extract_quadratic_coeffs") })
setGeneric("quad_form", function(object, expr) { standardGeneric("quad_form") })
## End of newly added generics
