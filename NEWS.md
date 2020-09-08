# CVXR 1.0-6-1.0-3

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

# CVXR pre-0.94-4

* Several wrong turns and much hand-wringing. Complete rewrite in
  preparation for release.
