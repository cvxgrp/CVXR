## Some crude tools for my use; will polish up later (BN)
## <CLASSNAME> class shadowing CPP class
<CLASSNAME> <- R6::R6Class("<CLASSNAME>",
                           private = list(
                               pkg = NA,
                               myClassName = NA,
                               ptr = NA
                           ),
                           active = list(
                               <ACTIVE_START>
                               <VAR> = function(value) {
                                   if (missing(value)) {
                                       rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                methodName = "get_<VAR>",
                                                                thisPkg = private$pkg)
                                       .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                                   } else {
                                       rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                                methodName = "set_<VAR>",
                                                                thisPkg = private$pkg)
                                       .Call(rcppFn, private$ptr, value, PACKAGE = private$pkg)
                                   }
                               }
                               <ACTIVE_END>
                           ),
                           public = list(
                               initialize = function(rcppClassName=gsub("CVXCanon.", class(self[1])) {
                                   private$pkg <- pkg <- getPackageName()
                                   private$myClassName <- rcppClassName
                                   private$ptr <- .Call(rcppMungedName(cppClassName = rcppClassName,
                                                                       methodName = "new",
                                                                       thisPkg = pkg),
                                                        PACKAGE = pkg)
                               }
                               <PUBLIC_START>
                               <NAME> = function() {
                                   rcppFn <- rcppMungedName(cppClassName = private$myClassName,
                                                            methodName = "<NAME>",
                                                            thisPkg = private$pkg)
                                   .Call(rcppFn, private$ptr, PACKAGE = private$pkg)
                               }
                               <PUBLIC_END>
                           ))
