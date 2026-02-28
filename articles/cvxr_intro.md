# Introduction to CVXR

## Overview

`CVXR` is an R package that provides an object-oriented modeling
language for convex optimization, similar to
[CVXPY](https://www.cvxpy.org/) for Python. It allows you to formulate
and solve convex optimization problems in a natural mathematical syntax.

## A Simple Example: Least Squares

Consider a simple linear regression problem where we want to estimate
parameters using a least squares criterion.

We generate synthetic data where we know the true model:

$$Y = X\beta + \epsilon$$

where $Y$ is a $100 \times 1$ vector, $X$ is a $100 \times 10$ matrix,
$\beta = ( - 4, - 3,\ldots,5)^{\top}$ is a $10 \times 1$ vector, and
$\epsilon \sim N(0,1)$.

``` r
set.seed(123)
n <- 100
p <- 10
beta <- -4:5

X <- matrix(rnorm(n * p), nrow = n)
Y <- X %*% beta + rnorm(n)
```

Using base R, we can estimate $\beta$ via `lm`:

``` r
ls.model <- lm(Y ~ 0 + X)
```

### The CVXR formulation

The same problem can be expressed as:

$$\underset{\beta}{\text{minimize}}\quad \parallel Y - X\beta \parallel_{2}^{2}$$

In CVXR, this translates directly:

``` r
library(CVXR)
#> 
#> Attaching package: 'CVXR'
#> The following objects are masked from 'package:stats':
#> 
#>     power, sd, var
#> The following objects are masked from 'package:base':
#> 
#>     norm, outer

betaHat <- Variable(p)
objective <- Minimize(sum((Y - X %*% betaHat)^2))
problem <- Problem(objective)
result <- psolve(problem)
```

The optimal value and estimated coefficients:

``` r
cat("Optimal value:", result, "\n")
#> Optimal value: 97.84759
cbind(CVXR = round(value(betaHat), 3),
      lm   = round(coef(ls.model), 3))
#>                lm
#> X1  -3.920 -3.920
#> X2  -3.012 -3.012
#> X3  -2.125 -2.125
#> X4  -0.867 -0.867
#> X5   0.091  0.091
#> X6   0.949  0.949
#> X7   2.076  2.076
#> X8   3.127  3.127
#> X9   3.961  3.961
#> X10  5.135  5.135
```

## Adding Constraints

The real power of CVXR is the ability to add constraints easily.

### Nonnegative Least Squares

Suppose we know the $\beta$s should be nonnegative:

``` r
problem <- Problem(objective, constraints = list(betaHat >= 0))
result <- psolve(problem)
round(value(betaHat), 3)
#>        [,1]
#>  [1,] 0.000
#>  [2,] 0.000
#>  [3,] 0.000
#>  [4,] 0.000
#>  [5,] 1.237
#>  [6,] 0.623
#>  [7,] 2.123
#>  [8,] 2.804
#>  [9,] 4.445
#> [10,] 5.207
```

### Custom Constraints

Now suppose $\beta_{2} + \beta_{3} \leq 0$ and all other $\beta$s are
nonnegative:

``` r
A <- matrix(c(0, 1, 1, rep(0, 7)), nrow = 1)
B <- diag(c(1, 0, 0, rep(1, 7)))

constraint1 <- A %*% betaHat <= 0
constraint2 <- B %*% betaHat >= 0

problem <- Problem(objective, constraints = list(constraint1, constraint2))
result <- psolve(problem, verbose = TRUE) ## verbose = TRUE for details
#> ─────────────────────────────── CVXR v1.8.0.9213 ───────────────────────────────
#> ℹ Problem: 1 variable, 2 constraints (QP)
#> ℹ Compilation: "OSQP" via CVXR::Dcp2Cone -> CVXR::CvxAttr2Constr -> CVXR::ConeMatrixStuffing -> CVXR::OSQP_QP_Solver
#> ℹ Compile time: 0.026s
#> ─────────────────────────────── Numerical solver ───────────────────────────────
#> -----------------------------------------------------------------
#>            OSQP v1.0.0  -  Operator Splitting QP Solver
#>               (c) The OSQP Developer Team
#> -----------------------------------------------------------------
#> problem:  variables n = 110, constraints m = 111
#>           nnz(P) + nnz(A) = 1210
#> settings: algebra = Built-in,
#>           OSQPInt = 4 bytes, OSQPFloat = 8 bytes,
#>           linear system solver = QDLDL v0.1.8,
#>           eps_abs = 1.0e-05, eps_rel = 1.0e-05,
#>           eps_prim_inf = 1.0e-04, eps_dual_inf = 1.0e-04,
#>           rho = 1.00e-01 (adaptive: 50 iterations),
#>           sigma = 1.00e-06, alpha = 1.60, max_iter = 10000
#>           check_termination: on (interval 25, duality gap: on),
#>           time_limit: 1.00e+10 sec,
#>           scaling: on (10 iterations), scaled_termination: off
#>           warm starting: on, polishing: on, 
#> iter   objective    prim res   dual res   gap        rel kkt    rho         time
#>    1   0.0000e+00   2.47e+01   4.44e+04  -6.89e+05   4.44e+04   1.00e-01    1.78e-04s
#>   50   1.1424e+02   3.47e+00   2.78e-04  -2.46e+02   3.47e+00   7.70e+00*   3.37e-04s
#>  125   1.2876e+03   6.48e-06   3.21e-06  -3.10e-03   6.48e-06   7.70e+00    5.88e-04s
#> plsh   1.2876e+03   5.82e-15   3.46e-13  -9.09e-13   3.46e-13   --------    6.96e-04s
#> 
#> status:               solved
#> solution polishing:   successful
#> number of iterations: 125
#> optimal objective:    1287.6297
#> dual objective:       1287.6297
#> duality gap:          -9.0949e-13
#> primal-dual integral: 6.8986e+05
#> run time:             6.96e-04s
#> optimal rho estimate: 1.41e+01
#> ──────────────────────────────────── Summary ───────────────────────────────────
#> ✔ Status: optimal
#> ✔ Optimal value: 1287.63
#> ℹ Compile time: 0.026s
#> ℹ Solver time: 0.001s
round(value(betaHat), 3)
#>         [,1]
#>  [1,]  0.000
#>  [2,] -2.845
#>  [3,] -1.711
#>  [4,]  0.000
#>  [5,]  0.664
#>  [6,]  1.178
#>  [7,]  2.329
#>  [8,]  2.414
#>  [9,]  4.212
#> [10,]  4.948
```

This demonstrates the chief advantage of CVXR: *flexibility*. Users can
quickly modify and re-solve a problem, making the package ideal for
prototyping new statistical methods. Its syntax is simple and
mathematically intuitive.

## Available Solvers

CVXR supports multiple solvers:

``` r
installed_solvers()
#> [1] "CLARABEL" "SCS"      "OSQP"     "HIGHS"
```

You can specify a solver explicitly:

``` r
psolve(problem, solver = "CLARABEL")
```

## Further Reading

- The [CVXR website](https://cvxr.rbind.io) has many worked examples
- The [CVXPY documentation](https://www.cvxpy.org/) covers the
  underlying mathematical framework
- Fu, Narasimhan, and Boyd (2020). “CVXR: An R Package for Disciplined
  Convex Optimization.” *Journal of Statistical Software*, 94(14).
  <doi:10.18637/jss.v094.i14>

## Session Info

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] CVXR_1.8.0.9213
#> 
#> loaded via a namespace (and not attached):
#>  [1] cli_3.6.5         knitr_1.51        rlang_1.1.7       xfun_0.56        
#>  [5] clarabel_0.11.2   highs_1.12.0-3    scs_3.2.7         textshaping_1.0.4
#>  [9] S7_0.2.1          jsonlite_2.0.0    backports_1.5.0   htmltools_0.5.9  
#> [13] gmp_0.7-5.1       ragg_1.5.0        sass_0.4.10       rmarkdown_2.30   
#> [17] grid_4.5.2        evaluate_1.0.5    jquerylib_0.1.4   fastmap_1.2.0    
#> [21] yaml_2.3.12       lifecycle_1.0.5   compiler_4.5.2    fs_1.6.6         
#> [25] Rcpp_1.1.1        osqp_1.0.0        systemfonts_1.3.1 lattice_0.22-7   
#> [29] digest_0.6.39     R6_2.6.1          checkmate_2.3.4   bslib_0.10.0     
#> [33] Matrix_1.7-4      tools_4.5.2       pkgdown_2.2.0     cachem_1.1.0     
#> [37] desc_1.4.3
```
