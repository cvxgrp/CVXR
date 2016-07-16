## CVXCanon.ProblemData class shadowing CPP class
CVXCanon.ProblemData <- R6::R6Class("CVXCanon.ProblemData",
                                    private = list(
                                        pkg = NA,
                                        myClassName = NA,
                                        ptr = NA
                                    ),
                                    active = list(
                                        V = function(value) {
                                            if (missing(value)) {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "get_V",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                            } else {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "set_V",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, value, PACKAGE = private$pkg)
                                            }
                                        }
                                       ,
                                        I = function(value) {
                                            if (missing(value)) {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "get_I",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                            } else {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "set_I",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, value, PACKAGE = private$pkg)
                                            }
                                        }
                                       ,
                                        J = function(value) {
                                            if (missing(value)) {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "get_J",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                            } else {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "set_J",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, value, PACKAGE = private$pkg)
                                            }
                                        }
                                       ,
                                        id_to_col = function(value) {
                                            if (missing(value)) {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "get_id_to_col",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                            } else {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "set_id_to_col",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, value, PACKAGE = private$pkg)
                                            }
                                        }
                                       ,
                                        const_to_row = function(value) {
                                            if (missing(value)) {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "get_const_to_row",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                            } else {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "set_const_to_row",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, value, PACKAGE = private$pkg)
                                            }
                                        }
                                       ,
                                        num_constraints = function(value) {
                                            if (missing(value)) {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "get_num_constraints",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                            } else {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "set_num_constraints",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, value, PACKAGE = private$pkg)
                                            }
                                        }
                                       ,
                                        vals = function(value) {
                                            if (missing(value)) {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "get_vals",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                            } else {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "set_vals",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, value, PACKAGE = private$pkg)
                                            }
                                        }
                                       ,
                                        row_idxs = function(value) {
                                            if (missing(value)) {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "get_row_idxs",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                            } else {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "set_row_idxs",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, value, PACKAGE = private$pkg)
                                            }
                                        }
                                       ,
                                        col_ptrs = function(value) {
                                            if (missing(value)) {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "get_col_ptrs",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                            } else {
                                                rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                         methodName = "set_col_ptrs",
                                                                         thisPkg = private$pkg)
                                                .Call(rcppFn, private$ptr, value, PACKAGE = private$pkg)
                                            }
                                        }
                                    ),
                                    public = list(
                                        initialize = function() {
                                            private$pkg <- pkg <- getPackageName()
                                            private$myClassName <- myClassName <- class(self)[1]
                                            private$ptr <- .Call(rcppMungedName(cppClassName = rcppClassName,
                                                                                methodName = "new",
                                                                                thisPkg = pkg),
                                                                 PACKAGE = pkg)
                                        }
                                       ,
                                        getV = function() {
                                            rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                     methodName = "getV",
                                                                     thisPkg = private$pkg)
                                            .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                        }
                                       ,
                                        getI = function() {
                                            rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                     methodName = "getI",
                                                                     thisPkg = private$pkg)
                                            .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                        }
                                       ,
                                        getJ = function() {
                                            rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                     methodName = "getJ",
                                                                     thisPkg = private$pkg)
                                            .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                        }
                                       ,
                                        getConstVec = function() {
                                            rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                     methodName = "getConstVec",
                                                                     thisPkg = private$pkg)
                                            .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                        }
                                    ))
