# CVXR 1.0-15

* Revert `clarabel` requirement and use enhance rather than import
  until really necessary ([Issue
  142](https://github.com/cvxgrp/CVXR/issues/142)).
* Update return codes for `user_limit` etc to be
  `infeasible_inaccurate` to match [`CVXPY` for gurobi](https://github.com/cvxpy/cvxpy/pull/1270)
* Add an S3 print method for result from `solve()`.
* Move `upper_tri_to_full` to C++
* Drop use of R6 classes

# CVXR 1.0-14

* Address further inefficiency in use of diagonal matrices in
  `qp2quad_form.R`
* Add initial interface to clarabel solver

# CVXR 1.0-13

* Address inefficient processing of cones for MOSEK (Issue 137
  reported by aszekMosek)
* Fix `extract_quadratic_coeffs` to use sparse matrix and sweep in
  place for better memory use (reported by Marissa Reitsma) 

# CVXR 1.0-12

* `Rmosek` to be removed from CRAN, so moved to drat repo
* Cleaned up some problematic Rd files shown in CRAN results

# CVXR 1.0-11

* Being more careful about coercing to `dgCMatrix`
  via `(as(as(<matrix>, "CsparseMatrix"), "generalMatrix"))`
* Modify all class inheritance checks to use `inherits()`

# CVXR 1.0-10

* Now requiring the updated scs 3.0 as an import

# CVXR 1.0-9

* Now importing ECOSolveR version 0.5.4 and higher
* Added fixes for Matrix 1.3 update
* Somewhat better typesetting for Power class documentation (Issue #86)

# CVXR 1.0-8 to 1.0-3

* Conforming to CRAN suggestions on non-CRAN packages, but no actual
  changes.

# CVXR 1.0-2

* Added exponential cone support for `MOSEK` and uncommented
  associated test.
* Added JSS publication DOI.

# CVXR 1.0-1

* Many small fixes to solver interfaces
* Reference semantics for `Parameter`
* Warm start for `OSQP` and updates to solver cache
* Solver parameter defaults explicit now
* New tests added for solver combinations

# CVXR 1.0

* Major release implementing reductions, many new solvers.

# CVXR 0.99

* Bug fix: duplicated integer and boolean indices. 
* Bug fix: correct typo in constraint specification to GLPK.
* Added tutorial articles based on v0.99 to [CVXR
website](https://cvxr.rbind.io) on using other solvers, integer
programming, MOSEK and GUROBI examples. 

# CVXR 0.98-1

* Minor typographical fixes.

# CVXR 0.98

* Dropped `delay_load` parameter dropped in
  `reticulate::import_from_path`, per changes in `reticulate`.

* Cleaned up hooks into reticulate for commercial solvers.

# CVXR 0.97-1

* Minor typo and documentation fixes.

# CVXR 0.97

* Added `LPSOLVE` via [`lpSolveAPI`](https://cran.r-project.org/package=lpSolveAPI)
* Added `GLPK` via [`Rglpk`](https://cran.r-project.org/package=Rglpk)
* Added `MOSEK` 
* Added `GUROBI`
* Bug fix: [issue #25](https://github.com/cvxgrp/CVXR/issues/25). 
  All CVXR expressions retain dimensions. Culprit was `drop =
  FALSE` (in function `Index.get_special_slice`) as suspected.  
  
# CVXR 0.96 

* Added a note that CVXR can probably be compiled from source for
  earlier versions of R. This is [issue
  #24](https://github.com/cvxgrp/CVXR/issues/24)

* Using [pkgdown](https://pkgdown.r-lib.org). This also addresses
  [issue #23](https://github.com/cvxgrp/CVXR/issues/23)

* Bug fix: [issue #28](https://github.com/cvxgrp/CVXR/issues/28)
  Function `intf_sign` (`interface.R`) was unnecessarily using a
  tolerance parameter, now eliminated.
  
# CVXR 0.95

* Updated Solver.solve to adapt to new ECOSolveR. Require version 0.4 of
  ECOSolveR now. 

* Updated `unpack_results` to behave exactly like in CVXPY. Added
  documentation and testthat tests. Documented in [Getting Faster
  Results article](https://cvxr.rbind.io/cvxr_examples/cvxr_speed/).
  

# CVXR 0.94-4

* First CRAN release 2017-11-20. 
