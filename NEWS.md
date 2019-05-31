# CVXR 0.99-6

* Mosek Glue modifications for MOSEK 8+, [issue 49](https://github.com/anqif/CVXR/issues/49).

* Format constraints bug for MOSEK (`constraints.R`) where `CVXOPT`
  and `MOSEK` were lumped together (thanks, Trevor Hastie).

# CVXR 0.99-5

* Bug fix for LogSumExp atom. This should address the issue reported
  on [StackOverflow](https://stackoverflow.com/questions/55737567/extension-to-the-cvxr-example-cvxr-kelly-strategy-not-dcp-compliant)
* Require ECOSolveR version 0.5.1 and above to avoid convolution
  example failure on 32-bit platforms.

# CVXR 0.99-4

* Updated tests to check against `scs` version 1.2.3, which now
  returns `optimal` rather than `optimal_inaccurrate` as it did
  earlier (`test-g01-non_optimal.R#35`).

* Some cleanup of manual tests. New SCS returns slightly different
  values for `Sigma` in manual test (`test-vignette`#646) for
  worst-case covariance example.

# CVXR 0.99-3

* Bug fix to MOSEK and GUROBI solver interfaces; dimensions now
  converted to expected python types.
* Sanitized src directory: proper naming of CVXR header to include
  while compiling attributes with Rcpp.

# CVXR 0.99-2

* Typo fixes to URLs.

# CVXR 0.99-1

* Bug fix: Updated python glue for both version 2 and 3
(`gurobiglue.py`)
* Bug fix: Workaround for zero extent sparse matrices in R not being
  handled by reticulate (`mosekglue.py`)

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
* Bug fix: [issue #25](https://github.com/anqif/CVXR/issues/25). 
  All CVXR expressions retain dimensions. Culprit was `drop =
  FALSE` (in function `Index.get_special_slice`) as suspected.  
  
# CVXR 0.96 

* Added a note that CVXR can probably be compiled from source for
  earlier versions of R. This is [issue
  #24](https://github.com/anqif/CVXR/issues/24)

* Using [pkgdown](https://pkgdown.r-lib.org). This also addresses
  [issue #23](https://github.com/anqif/CVXR/issues/23)

* Bug fix: [issue #28](https://github.com/anqif/CVXR/issues/28)
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
