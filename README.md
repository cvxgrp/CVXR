# CVXR

An R modeling language for convex optimization problems.

## Note

This package is still under active development. You can install it, but crucial glue to the Canon interface which bridges any language---in our case R---and C-level libraries still requires work to do anything useful.

## Installation

To install, make sure you have the development tools for R available, including the C compilers etc.

1. Install the packages `devtools` `Rcpp`, `RcppEigen`, `uuid`, `bitops`, `hashmap`, `rstackdeque`,
2. In R, execute
```
library(devtools)
install("anqif/cvxr")
```

## Support

You may experiment with the DSL expressing problems using the current version.
We are unable to provide any support until the package is officially released.




