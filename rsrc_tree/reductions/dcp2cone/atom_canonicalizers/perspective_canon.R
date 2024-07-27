## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/perspective_canon.py

Dcp2Cone.perspective_canon <- function(expr, args) {
  # Only working for minimization right now.

  if(is_convex(expr@f))
    aux_prob <- Problem(Minimize(expr@f))
  else
    aux_prob <- Problem(Maximize(expr@f))
  # Does numerical solution value of epigraph t coincide with expr.f numerical
  # value at opt?
  solver_opts <- list("use_quad_obj" = FALSE)
  chain <- .construct_chain(aux_prob, solver_opts=solver_opts, ignore_dpp=TRUE)
  num_red <- length(chain@reductions)
  if(num_red > 1)
    chain@reductions = chain@reductions[[seq(num_red - 1)]]  # skip solver reduction
  prob_canon <- apply(chain, aux_prob)[[1]]  # grab problem instance
  # get cone representation of c, A, and b for some problem.

  c <- as.vector(prob_canon@c)
  c_len <- length(c)
  c <- c[seq_len(c_len - 1)]
  d <- as.vector(prob_canon@c)
  d <- d[c_len]

  Ab <- matrix(prob_canon@A, ncol = length(c) + 1, byrow = FALSE)
  Ab_ncol <- ncol(Ab)
  A <- Ab[, seq_len(Ab_ncol - 1)]
  b <- Ab[, Ab_ncol]

  # given f in epigraph form, aka epi f = \{(x,t) | f(x) \leq t\}
  # = \{(x,t) | Fx +tg + e \in K} for K a cone, the epigraph of the
  # perspective, \{(x,s,t) | sf(x/s) \leq t} = \{(x,s,t) | Fx + tg + se \in K\}
  # If I have the problem "minimize f(x)" written in the CVXPY compatible
  # "c^Tx, Ax+b \in K" form, I can re-write this in the graph form above via
  # x,t \in \epi f iff Ax + b \in K and t-c^Tx \in R_+ which I can further write
  # with block matrices as Fx + tg + e \in K \times R_+
  # with F = [A ], g = [0], e = [b]
  #          [-c]      [1]      [-d]

  # Actually, all we need is Ax + 0*t + sb \in K, -c^Tx + t - ds >= 0

  t <- Variable()
  s <- args[[1]]
  x_canon <- prob_canon@x
  constraints <- list()

  if(!is.null(dim(A))) {
    # Rules out the case where f is affine and requires no additional
    # constraints.
    x_pers <- A@x_canon + s*b

    i <- 1
    for(con in prob_canon@constraints) {
      sz <- size(con)
      var_slice <- x_pers[seq(i, i + sz - 1)]
      pers_constraint <- form_cone_constraint(var_slice, con)
      constraints <- c(constraints, list(pers_constraint))
      i <- i + sz
    }
  }

  constraints <- c(constraints, list(-c@x_canon + t - s*d >= 0))

  # recover initial variables

  end_inds <- c(sort(unlist(prob_canon@var_id_to_col), decreasing = FALSE), nrow(x_canon) + 1)

  for(var in variables(expr@f)) {
    start_ind <- prob_canon@var_id_to_col[[as.character(id(var))]]
    end_ind <- end_inds[which(end_inds == start_ind) + 1]

    if(var@attributes$diag)  # checking for diagonal first because diagonal is also symmetric
      constraints <- c(constraints, list(diag(var) == x_canon[start_ind:(end_ind - 1)]))
    else if(is_symmetric(var) && size(var) > 1) {
      n <- nrow(var)
      inds <- which(upper.tri(matrix(1, nrow = n, ncol = n), diag = TRUE), arr.ind = TRUE)  # includes diagonal
      constraints <- c(constraints, list(var[inds] == x_canon[start_ind:(end_ind - 1)]))
    } else
      constraints <- c(constraints, list(vec(var) == x_canon[start_ind:(end_ind - 1)]))
  }

  if(is_convex(expr@f))
    return(list(t, constraints))
  else
    return(list(-t, constraints))
}
