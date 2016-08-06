## Solution class shadowing CPP class
CVXcanon <- R6::R6Class("CVXcanon.CVXcanon",
                        private = list(
                            pkg = NA,
                            myClassName = NA,
                            solution = NA
                        )
                       ,
                        active = list(
                        )
                       ,
                        public = list(
                            initialize = function() {
                                private$pkg <- getPackageName()
                                private$myClassName <- class(self)[1]
                            }
                           ,
                            solve = function(sense, objective, constraints, solverOptions) {
                                pkg <- private$pkg
                                ptr <- .Call(rcppMungedName(cppClassName = private$myClassName,
                                                            methodName = "solve",
                                                            thisPkg = pkg),
                                             sense, objective$getXPtr(), constraints$getXPtr(),
                                             solverOptions,
                                             PACKAGE = pkg)
                                CVXcanon.Solution$new(ptr)
                            }
                        ))
