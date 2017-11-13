## An R6 List class with reference semantics mainly for keeping objects in scope
#' @importFrom R6 R6Class
R6List <- R6::R6Class("R6List",
                      private = list(
                          r6list = NA
                      )
                     ,
                      public = list(
                          initialize = function() {
                              private$r6list <- list()
                          }
                         ,
                          getList = function() {
                              private$r6list
                          }
                          ,
                          append = function(what) {
                              n <- length(private$r6list)
                              private$r6list[[n+1]] <- what
                          }
                         ,
                          toString = function() {
                              result <- sapply(private$r6list, function(x) x$toString())
                              result <- paste(result, collapse = ", ")
                              sprintf("[ %s ]", result)
                          }
                          ,
                          print = function() {
                              print(self$toString())
                          }
                      ))
