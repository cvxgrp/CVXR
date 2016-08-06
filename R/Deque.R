Deque <- R6::R6Class("Deque",
                    private = list(
                        queue = NA
                    ),
                    public = list(
                        initialize = function(init = list()) {
                            private$queue <- init
                        },
                        ## popleft, remove and return an element from the left side of the queue,
                        ## raise error if no element.
                        popleft = function() {
                            if (length(private$queue) > 0) {
                                result <- private$queue[[1]]
                                private$queue <- private$queue[-1]
                                result
                            } else {
                                stop("Deque: empty queue")
                            }
                        },
                        ## append, adds to the right end of the queue
                        ##
                        append = function(what) {
                            n <- length(private$queue)
                            private$queue[[n+1]] <- what
                            invisible(private$queue)
                        },
                        length = function() {
                            length(private$queue)
                        },
                        getQ = function() private$queue,
                        print = function() {
                            print(private$queue)
                        }
                    ))

