# CVXR 1.8.1 — CRAN Submission Comments

## Overview

This is a major rewrite of CVXR, migrating from S4 to S7 and aligning
the codebase with CVXPY 1.8. The rewrite delivers ~4-5x faster
performance, improved maintainability, and full feature parity with
CVXPY including DPP, DGP, DQCP, complex variables, and 13 solver
backends.

## R CMD check results

0 errors | 0 warnings | 0 notes

Tested on:
- macOS Tahoe 26.3 (aarch64), R 4.5.2
- (additional platforms via GitHub Actions)

## Reverse dependency check

Checked all 37 reverse dependencies (36 CRAN + 1 Bioconductor).

**Results**: 5 OK, 7 install failures, 19 broken, 6 warning-only.

All failures stem from intentional API changes in this major version:

| Breaking change | Packages affected | Severity |
|:----------------|:-----------------:|:---------|
| `CVXR::solve` / `importFrom(CVXR, solve)` no longer available | 24 | 7 install failures, 10 runtime errors, 7 warnings |
| `Variable(rows, cols)` → `Variable(c(rows, cols))` | 3 | Runtime errors |
| `CVXR::diag` no longer exported (use `diag()` or `DiagVec()`/`DiagMat()`) | 2 | Runtime errors / warnings |
| `CVXR::diff` no longer exported (use `diff()` or `cvxr_diff()`) | 2 | Warnings |
| `solve()` return value format change | 3 | Runtime errors |
| `sd`/`var`/`outer` export masking | 3 | Warnings |
| `getValue()` deprecation | 1 | Warning |
| S7 object coercion | 2 | Runtime errors |

### Migration path

Migration path was specified in email to package maintainers as well as the [CVXR website](https://cvxr.rbind.io/whatsnew#migration-guide). 

The primary change is `solve(prob)` → `psolve(prob)`. A backward-compatible
`solve()` S3 method is registered at runtime (works for attached usage),
but `CVXR::solve()` and `importFrom(CVXR, solve)` no longer resolve
because `solve` is not an explicit NAMESPACE export.

All affected maintainers have been notified with specific, per-package
migration instructions. The changes required are minimal (typically 1-3
lines per package).

### Unaffected packages (5)

DebiasInfer, filling, mlr3fairness, portfolioBacktest, SIHR

### Affected packages by issue

**Install failures** (7): ANCOMBC, fungible, migest, PlackettLuce,
Riemann, spBPS, wdnet

**Runtime errors** (12): aramappings, ccar3, cuadramelo, DiSCos,
EmpiricalDynamics, EpiForsk, glmmrOptim, hedgedrf, kantorovich,
pvEBayes, rclsp, scpi, SLSEdesign, spStack, starnet, tramnet, transreg

**Warning only** (6): fairml, graphicalExtremes, HonestDiD,
MaximinInfer, PortfolioAnalytics, Rdimtools, RobustIV, WRI
