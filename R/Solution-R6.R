## Solution class shadowing CPP class
CVXCanon.Solution <- R6::R6Class("CVXCanon.Solution",
                                 private = list(
                                     pkg = NA,
                                     myClassName = NA,
                                     ptr = NA
                                 ),
                                 active = list(
                                     status = function(value) {
                                         if (missing(value)) {
                                             rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                      methodName = "get_status",
                                                                      thisPkg = private$pkg)
                                             .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                         } else {
                                             rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                      methodName = "set_status",
                                                                      thisPkg = private$pkg)
                                             .Call(rcppFn, private$ptr, value, PACKAGE = private$pkg)
                                         }
                                     }
                                    ,
                                     optimal_value = function(value) {
                                         if (missing(value)) {
                                             rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                      methodName = "get_optimal_value",
                                                                      thisPkg = private$pkg)
                                             .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                         } else {
                                             rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                      methodName = "set_optimal_value",
                                                                      thisPkg = private$pkg)
                                             .Call(rcppFn, private$ptr, value, PACKAGE = private$pkg)
                                         }
                                     }
                                    ,
                                     primal_values = function(value) {
                                         if (missing(value)) {
                                             rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                      methodName = "get_primal_values",
                                                                      thisPkg = private$pkg)
                                             .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                         } else {
                                             rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                      methodName = "set_primal_values",
                                                                      thisPkg = private$pkg)
                                             .Call(rcppFn, private$ptr, value, PACKAGE = private$pkg)
                                         }
                                     }
                                    ,
                                     dual_values = function(value) {
                                         if (missing(value)) {
                                             rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                      methodName = "get_dual_values",
                                                                      thisPkg = private$pkg)
                                             .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                         } else {
                                             rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                      methodName = "set_dual_values",
                                                                      thisPkg = private$pkg)
                                             .Call(rcppFn, private$ptr, value, PACKAGE = private$pkg)
                                         }
                                     }
                                 ),
                                 public = list(
                                     initialize = function() {
                                         private$pkg <- pkg <- getPackageName()
                                         private$myClassName <- myClassName <- class(self)[1]
                                         private$ptr <- .Call(rcppMungedName(cppClassName = myClassName,
                                                                             methodName = "new",
                                                                             thisPkg = pkg),
                                                              PACKAGE = pkg)
                                     }
                                 ))
