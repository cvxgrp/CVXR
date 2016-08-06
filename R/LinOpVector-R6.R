## CVXcanon.LinOpVector class shadowing CPP class
CVXcanon.LinOpVector <- R6::R6Class("CVXcanon.LinOpVector",
                                    private = list(
                                        pkg = NA,
                                        myClassName = NA,
                                        linOps = NA,
                                        ptr = NA ## the rcpp XPtr
                                    ),
                                    active = list(
                                    ),
                                    public = list(
                                        initialize = function() {
                                            private$linOps <- list()
                                            private$pkg <- pkg <- getPackageName()
                                            private$myClassName <- myClassName <- class(self)[1]
                                            private$ptr <- .Call(rcppMungedName(cppClassName = myClassName,
                                                                                methodName = "new",
                                                                                thisPkg = pkg),
                                                                 PACKAGE = pkg)
                                        }
                                       ,
                                        getXPtr = function() {
                                            private$ptr
                                        }
                                       ,
                                        push_back = function(R6LinOp) {
                                            n <- length(private$linOps)
                                            private$linOps[[n+1]] <- R6LinOp
                                            ## Needs modification by hand for arguments
                                            rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                     methodName = "push_back",
                                                                     thisPkg = private$pkg)
                                            .Call(rcppFn, private$ptr , R6LinOp$getXPtr(), PACKAGE = private$pkg)
                                        }
                                       ,
                                        toString = function() {
                                            result <- sapply(private$linOps, function(x) x$toString())
                                            result <- paste(result, collapse = ", ")
                                            sprintf("[ %s ]", result)
                                        }
                                        ,
                                        print = function() {
                                            print(self$toString())
                                        }

                                    ))
