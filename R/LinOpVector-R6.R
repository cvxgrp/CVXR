## LinOpVector class shadowing CPP class
CVXcanon.LinOpVector <- R6::R6Class("CVXcanon.LinOpVector",
                                    private = list(
                                        linOps = NA,
                                        ptr = NA ## the rcpp XPtr
                                    ),
                                    active = list(
                                    ),
                                    public = list(
                                        initialize = function() {
                                            private$linOps <- list()
                                            private$ptr <- .Call("_CVXR_LinOpVector__new", PACKAGE = "CVXR")
                                        }
                                       ,
                                        getXPtr = function() {
                                            private$ptr
                                        }
                                       ,
                                        getList = function() {
                                            private$linOps
                                        }
                                       ,
                                        push_back = function(R6LinOp) {
                                            n <- length(private$linOps)
                                            private$linOps[[n+1]] <- R6LinOp
                                            ## Needs modification by hand for arguments
                                            .Call("_CVXR_LinOpVector__push_back", private$ptr , R6LinOp$getXPtr(), PACKAGE = "CVXR")
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
