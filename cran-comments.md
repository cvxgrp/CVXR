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
- macOS Tahoe 26.3.1 (aarch64), R 4.5.2
- (additional platforms via GitHub Actions)

## Reverse dependency check

**Date**: 2026-03-05

Checked all 37 reverse dependencies (36 CRAN + 1 Bioconductor).

**Results**: 5 OK, 7 install failures, 25 with new problems (17
runtime errors, 8 warning-only).

All failures stem from intentional API changes in this major version:

| Breaking change | Packages affected | Severity |
|:----------------|:-----------------:|:---------|
| `solve` no longer exported (`importFrom(CVXR, solve)` / `CVXR::solve()` fail) | 28 | 7 install failures, 14 runtime errors, 7 warnings |
| `Variable(rows, cols)` → `Variable(c(rows, cols))` | 5 | Runtime errors (aramappings, ccar3, cuadramelo, Riemann, wdnet) |
| `psolve()` returns scalar (not result list) | 2 | Runtime errors (kantorovich, spStack) |
| `CVXR::diag` no longer exported (use bare `diag()` or `DiagVec()`/`DiagMat()`) | 2 | Runtime errors / warnings |
| `CVXR::diff` no longer exported (use bare `diff()` or `cvxr_diff()`) | 2 | Warnings |
| `sd`/`var`/`outer` export masking | 3 | Warnings |
| Matrix S4 `%*%` CVXR S7 dispatch conflict (fix: `as_cvxr_expr()`) | 6 | 3 confirmed crashes (rclsp, tramnet, glmmrOptim) + 3 conditional (scpi, PlackettLuce, PortfolioAnalytics) |
| `sum_entries(x, axis=2)` shape `(1,n)` vs old 1D (needs `t()` on target) | 1 | Runtime infeasibility (wdnet) |
| `getValue()` deprecation | 1 | Warning |

### Migration path

The primary change is `solve(prob)` → `psolve(prob)` + `value(var)`.
A backward-compatible `solve()` S3 method is registered at runtime
(works for attached usage), but `CVXR::solve()` and
`importFrom(CVXR, solve)` no longer resolve because `solve` is not an
explicit NAMESPACE export.

### Fix verification

6 of 7 install-failure packages were patched and verified (PlackettLuce
excluded — needs additional Matrix interop fix):

| Package | Install | Tests | Notes |
|:--------|:-------:|:-----:|:------|
| ANCOMBC | PASS | N/A | Examples need Bioconductor Suggests |
| fungible | PASS | PASS | |
| migest | PASS | PASS | |
| Riemann | PASS | PASS | Also needed `Variable(m,n)` → `Variable(c(m,n))` |
| spBPS | PASS | PASS | |
| wdnet | PASS | 74/74 | Also needed Variable sig + `t()` on `sum_entries` axis=2 target |

All 30 affected maintainers have been notified with specific,
per-package migration instructions including exact file:line
references and code diffs. The changes required are minimal (typically
1-5 lines per package). Per-package details are in
`notes/revdep_messages/`.

### Unaffected packages (5)

DebiasInfer, filling, mlr3fairness, portfolioBacktest, SIHR

### Affected packages by issue

**Install failures** (7 — all `importFrom(CVXR, solve)`):
ANCOMBC, fungible, migest, PlackettLuce (also Matrix S4 %*% risk),
Riemann (also Variable sig), spBPS, wdnet (also Variable sig +
axis shape)

**Runtime errors** (17):
aramappings (Variable sig + solve), ccar3 (Variable sig + solve),
cuadramelo (Variable sig + solve), DiSCos (solve), EmpiricalDynamics
(solve → getValue), EpiForsk (solve), glmmrOptim (solve + Matrix S4
%*% S7), hedgedrf (solve), kantorovich (psolve result API), pvEBayes
(solve), rclsp (Matrix S4 %*% S7 + solve API), scpi (solve + sd/var
masking + Matrix S4 %*% risk), SLSEdesign (diag + solve), spStack
(psolve result API), starnet (solve), tramnet (solve falls to
base::solve + Matrix S4 %*% S7), transreg (diff + solve)

**Warning only** (8):
fairml (solve), graphicalExtremes (diag + solve), HonestDiD (solve),
MaximinInfer (getValue deprecation), PortfolioAnalytics (solve +
Matrix S4 %*% risk), Rdimtools (sd/var masking), RobustIV (solve),
WRI (diff)
