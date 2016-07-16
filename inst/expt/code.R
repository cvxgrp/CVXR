Deque <- R6Class("Deque",
                 private = list(
                     queue = NA
                 ),
                 public = list(
                     initialize = function() {
                         private$queue = list()
                     },
                     pop = function() {
                         if (length(private$queue) > 0) {
                             result <- private$queue[[1]]
                             private$queue <- private$queue[-1]
                             result
                         } else {
                             stop("Deque: empty queue")
                         }
                     },
                     push_back = function(what) {
                         n <- length(private$queue)
                         private$queue[[n+1]] <- what
                     },
                     length = function() {
                         length(private$queue)
                     },
                     getQ = function() private$queue,
                     print = function() {
                         print(private$queue)
                     }
                 ))



build_lin_op_tree <- function(root_linR, tmp, verbose = FALSE) {
    Q <- Deque$new()
    root_linC <- CVXCanon.LinOp$new()
    Q$push_back(list(linR = root_linR, linC = root_linC))

    while(Q$length() > 0) {
        node <- Q$pop()
        linR <- node$linR
        linC <- node$linC

        ## Updating the arguments our LinOp
        ## tmp is a list
        for(argR in linR$args) {
            tree <- CVXCanon.LinOp$new()
            tmp <- c(tmp, tree)
            Q$push_back(list(linR = argR, linC = tree))
            linC$args_push_back(tree)
        }

        ## Setting the type of our LinOp; at the C level, it is an ENUM!
        ## Can we avoid this case conversion and use UPPER CASE to match C?
        linC$type <- linop_type2Int[[toupper(linR$type)]] ## Check with Anqi

        ## Setting size
        linC$size_push_back(as.integer(linR$size[1]))
        linC$size_push_back(as.integer(linR$size[2]))

        ## Loading the problem data into the approriate array format
        if(!is.na(linR$data)) {
            if(length(linR$data) == 2 && is(linR$data[1], 'slice'))
                ## ASK Anqi about this
                set_slice_data(linC, linR)  ## TODO
            else if(is.numeric(linR$data) || is.integer(linR$data))
                linC$dense_data <- format_matrix(linR$data, 'scalar')
            else if(linR$data$class == 'LinOp' && linR$data$type == 'scalar_const')
                linC$dense_data <- format_matrix(linR$data$data, 'scalar')
            else
                set_matrix_data(linC, linR)
        } else
            set_matrix_data(linC, linR)
    }

    root_linC
}

