## An R6 class to hold references to objects so that they are not garbage collected until
## this class is destroyed
## imports UUIDgenerate::UUIDG
ReferenceHolder <- R6::R6Class("ReferenceHolder",
                               private = list(
                                   env = NA
                               )
                              ,
                               public = list(
                                   initialize = function() {
                                       private$env <- = new.env(parent = emptyenv())
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
                          print = function() {
                              for (x in private$r6list) {
                                  print(x)
                              }
                          }
                      ))
