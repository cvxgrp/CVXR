# CVXR 1.9.2 — CRAN Submission Comments

## Overview

This is a **correctness release**. It tracks CVXPY 1.9.2 and fixes several
bugs that returned plausible-looking wrong answers with no error and no
warning. The most consequential has been present since 1.8.2:

* **A constant weight on `sum(square(x))` was discarded when building the
  quadratic objective**, so penalized least-squares problems were solved with
  the wrong penalty. Ridge regression is the common case — every value of
  `lambda` returned the same fit. The solution returned was optimal for a
  different problem than the one written. `sum_squares()`, `quad_form()` and
  `p_norm(x, 2)^2` were unaffected, which is why it went unnoticed: a pure
  quadratic has the same minimizer under any positive weight, so the defect is
  invisible unless the objective also has a linear term.
* **A rejected solve left a stale compiled chain on the problem object**, so
  asking again could silently return an answer computed a different way (for
  example, `gp = TRUE` on a non-DGP problem erroring once and then solving as
  an ordinary DCP problem).
* **`sqrt()` and `inv_pos()` canonicalized to the wrong cone**, emitting
  `PowCone3D` where CVXPY emits an SOC. SOC-only solvers rejected such problems
  outright, and one solver returned a silently wrong answer on the DQCP form.
* A complex `quad_form()` with a negative semidefinite matrix passed `is_dcp()`
  and then failed to compile.

The release also completes a **member-level completeness audit** of CVXR
against CVXPY 1.9.2 (nine findings, all fixed), and is **about a third faster
than 1.9.1** with substantially lower memory use on larger problems.

## R CMD check results

**0 errors | 0 warnings | 0 notes**

Tested on:
- macOS 26.6.2 (aarch64-apple-darwin23), R 4.6.1 (2026-06-24) — `--as-cran`,
  64 checks OK
- Test suite: 16658 passing, 0 failing

### Sanitizers (UBSAN + ASAN)

The package contains compiled C++, and a macOS check does not exercise the
sanitizer classes (nor the pragma scan). It was therefore additionally checked
on **Linux R-devel built with ASAN + UBSAN** (`rocker/r-devel-san`,
R-devel 2026-08-17 r90424, x86_64), against the 2026-08-17 CRAN-form tarball
of this release. The only compiled-code change since that run is a three-line
NaN canonicalization in a constant-keying helper (`src/RcppCseKeys.cpp`,
2026-08-23): the helper itself runs on every solve and was fully exercised
under the sanitizers; the new branch adds only `ISNAN`/`R_IsNA` tests and a
constant store, no new memory or arithmetic operations. Results of that run:

* package **installs and loads** under `-fsanitize=undefined,bounds-strict`
* `checking tests ... OK` — 1139 assertions, 0 failures
* **0 `runtime error:` lines and 0 AddressSanitizer reports**

The findings on that platform are all container artifacts rather than package
defects: the dependency WARNING is the absent optional solvers (see below),
and the NOTEs are a missing `pandoc`, a missing HTML `tidy`, no vignettes in
the CRAN build, and one example (`psolve`, 5.2s) that exceeds 5s only because
sanitizer instrumentation and CPU emulation slow it by roughly two orders of
magnitude; it runs in well under 0.1s natively.

Scope, stated precisely: this exercises the CRAN test subset, so it covers
CVXR's own compiled code on every path that subset reaches. The commercial
solver interfaces (MOSEK, Gurobi, CPLEX, XPRESS) cannot be installed on that
platform and were not exercised.

## User-visible changes

These are deliberate and documented in NEWS. None is expected to affect a
reverse dependency, and the reverse-dependency results below bear that out.

* `dual_value()` now returns a matrix shaped like its constraint rather than a
  flat vector, which makes its documented return type accurate — it previously
  promised a matrix and returned a vector. Linear indexing of a dual (`d[3]`)
  is unaffected, since R matrices index linearly in column-major order; code
  calling `length()` or `dim()` on a dual may need adjusting. Scalar
  constraints, `SOC` and `ExpCone` are deliberately unchanged, matching CVXPY.
* A parameterized problem that is not DPP now emits a warning saying that
  subsequent solves will not be faster, matching CVXPY. This is a new warning
  on problems that previously solved silently; the answers are unchanged.
* Mixed-integer problems using a *partial* `boolean`/`integer` index list now
  correctly require a mixed-integer solver. A problem that previously appeared
  to solve with a continuous solver now reports that no suitable solver is
  available if none is installed.
* Five error conditions (`DPPError`, `DGPError`, `DQCPError`, `DNLPError`,
  `ParameterError`) are now catchable by class, as in CVXPY. Previously they
  could only be matched by message text.

## Reverse dependency check

**Date**: 2026-08-23 · `revdepcheck` 1.0.0.9002, R 4.6.1, macOS 26.6.2
(aarch64), Bioconductor included. This is a fresh run against the final
submission tarball (an earlier full run on 2026-08-20 gave the same result).

We checked **all 39 reverse dependencies**, comparing `R CMD check` results
across CVXR 1.9.1 (CRAN) and 1.9.2 (this submission).

* **0 packages broken**
* **0 packages failed to check**

`revdepcheck` initially reported one new problem, `aramappings`. It is a
**false positive from the parallel harness, not a CVXR regression**, and it was
verified as such rather than assumed:

```
> cl <- parallel::makeCluster(NCORES)
Error in serverSocket(port = port) :
  creation of server socket failed: port 11219 cannot be opened
```

The check ran with 15 concurrent workers, and `aramappings`' examples start
their own PSOCK cluster; the combination exhausts the local port space. The
failure is in `parallel::makeCluster`, before any CVXR code is reached — and
its "newly fixed" column shows the same socket error against CRAN CVXR 1.9.1,
confirming the artifact is bidirectional harness noise. Re-checked
**serially** against the same CVXR 1.9.2 library that revdepcheck built, the
package is `Status: OK`. (The same package failed the same way in the
2026-08-12 and 2026-08-20 runs at the same worker count, and was likewise
clean serially each time.)

No reverse dependency required a change for this release. This is consistent
with the nature of the release: the user-visible changes listed above either
correct wrong answers or add a warning, and none alters a function signature.

## Previous submission

CVXR 1.9.1 is the current CRAN version. The 1.8.1 submission was the S4-to-S7
rewrite and carried extensive API-break notes; those breaks are long settled
and no longer apply. This release makes no such changes.
