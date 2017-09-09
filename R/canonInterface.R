## Added format_matrix and set_matrix_data.
get_problem_matrix <- function(constrs, id_to_col = integer(0), constr_offsets = integer(0)) {
    cvxCanon <- CVXcanon$new()
    linOps <- lapply(constrs, function(constr) { constr$expr })
    lin_vec <- CVXcanon.LinOpVector$new()

    ## KLUDGE: Anqi, fix id_to_col to have proper names!
    # if (is.null(names(id_to_col))) names(id_to_col) <- unlist(id_to_col)
    ## END OF KLUDGE

    id_to_col_C <- id_to_col
    ## The following ensures that id_to_col_C is an integer vector
    ## with names retained. This is the C equivalent of map<int, int> in R
    storage.mode(id_to_col_C) <- "integer"

    ##if (any(is.na(id_to_col)))
    ##    id_to_col <- c()

    ## Loading the variable offsets from our R list into a C++ map

    ## for (id in names(id_to_col)) {
    ##     col <- id_to_col[[id]]
    ##     id_to_col_C$map(key = as.integer(id), value = as.integer(col))
    ## }

    ## This array keeps variables data in scope after build_lin_op_tree returns
    tmp <- R6List$new()
    for (lin in linOps) {
        tree <- build_lin_op_tree(lin, tmp)
        tmp$append(tree)
        lin_vec$push_back(tree)
    }

    ## REMOVE this later when we are sure
    if (typeof(constr_offsets) != "integer") {
        stop("get_problem_matrix: expecting integer vector for constr_offsets")
    }

    if (length(constr_offsets) == 0)
        problemData <- cvxCanon$build_matrix(lin_vec, id_to_col_C)
    else {
        ## Load constraint offsets into a C++ vector
        ##constr_offsets_C <- CVXCanon.IntVector$new()
        ##for (offset in constr_offsets)
        ##    constr_offsets_C$push_back(as.integer(offset))
        constr_offsets_C <- constr_offsets
        storage.mode(constr_offsets_C) <- "integer"
        problemData <- cvxCanon$build_matrix(lin_vec, id_to_col_C, constr_offsets_C)
    }

    ## Unpacking
    ## V <- problemData$getV()
    ## I <- problemData$getI()
    ## J <- problemData$getJ()
    ## const_vec <- problemData$getConstVec()

    list(V = problemData$getV(), I = problemData$getI(), J = problemData$getJ(),
         const_vec = matrix(problemData$getConstVec(), ncol = 1))
}

format_matrix <- function(matrix, format='dense') {
    ## Returns the matrix in the appropriate form,
    ## so that it can be efficiently loaded with our swig wrapper
    if (format == 'dense') {
        ## Ensure is 2D.
        as.matrix(matrix)
    } else if (format == 'sparse') {
        Matrix::Matrix(matrix, sparse = TRUE)
    } else if (format == 'scalar') {
        ## Should this be a 1x1 matrix?  YESSSSS as I later found out.
        as.matrix(matrix)
    } else {
        stop(sprintf("format_matrix: format %s unknown", format))
    }
}


set_matrix_data <- function(linC, linR) {
    ## Calls the appropriate CVXcanon function to set the matrix
    ## data field of our C++ linOp.

    if (is.list(linR$data) && linR$data$class == "LinOp") {
        if (linR$data$type == 'sparse_const') {
            linC$sparse_data <- format_matrix(linR$data$data, 'sparse')
        } else if (linR$data$type == 'dense_const') {
            linC$dense_data <- format_matrix(linR$data$data)
        } else {
            stop(sprintf("set_matrix_data: data.type %s unknown", linR$data$type))
        }
    } else {
        if (linR$type == 'sparse_const') {
            linC$sparse_data <- format_matrix(linR$data, 'sparse')
        } else {
            linC$dense_data <- format_matrix(linR$data)
        }
    }
}

set_slice_data <- function(linC, linR) {  ## What does this do?
    for (i in seq.int(length(linR$data) - 1L)) {  ## the last element is "class"
        sl <- linR$data[[i]]
        ## In R this is just a vector of ints

        ## ## Using zero based indexing throughout
        ## start_idx <- 0L
        ## if (!is.na(sl$start_idx))
        ##     start_idx <- sl$start_idx - 1L  ## Using zero-based indexing

        ## stop_idx <- linR$args[[1]]$size[i] - 1L
        ## if (!is.na(sl$stop_idx))
        ##     stop_idx <- sl$stop_idx - 1L

        ## step <- 1L
        ## if(!is.na(sl$step))
        ##     step <- sl$step

        ## ## handle [::-1] case
        ## if(step < 0 && is.na(sl$start_idx) && is.na(sl$stop_idx)) {
        ##     tmp <- start
        ##     start_idx <- stop_idx - 1
        ##     stop_idx <- tmp
        ## }

        ##for(var in c(start_idx, stop_idx, step))
        ##    vec$push_back(var)
        ## vec <- c(start_idx, stop_idx, step)
        if (length(sl) == 1L) {
            vec <- c(sl - 1L, sl, 1L)
        } else if (length(sl) == 2L) {
            vec <- c(sl[1L] - 1L, sl[2L], 1L)  # Using zero-based indexing, and step assumed to be 1.
        } else {
            r <- range(sl)
            vec <- c(r[1L] - 1L, r[2L], 1L)
        }

        ##vec <- c(sl, 1L)  # Using 1-based indexing, and step assumed to be 1.
        linC$slice_push_back(vec)
    }
}

build_lin_op_tree <- function(root_linR, tmp, verbose = FALSE) {
    Q <- Deque$new()
    root_linC <- CVXcanon.LinOp$new()
    Q$append(list(linR = root_linR, linC = root_linC))

    while(Q$length() > 0) {
        node <- Q$popleft()
        linR <- node$linR
        linC <- node$linC

        ## Updating the arguments our LinOp
        ## tmp is a list
        for(argR in linR$args) {
            tree <- CVXcanon.LinOp$new()
            tmp$append(tree)
            Q$append(list(linR = argR, linC = tree))
            linC$args_push_back(tree)
        }

        ## Setting the type of our LinOp; at the C level, it is an ENUM!
        ## Can we avoid this case conversion and use UPPER CASE to match C?
        linC$type <- toupper(linR$type) ## Check with Anqi

        ## Setting size
        linC$size_push_back(as.integer(linR$size[1]))
        linC$size_push_back(as.integer(linR$size[2]))

        ## Loading the problem data into the approriate array format
        if(!is.null(linR$data)) {
            ## if(length(linR$data) == 2 && is(linR$data[1], 'slice'))
            if (length(linR$data) == 3L && linR$data[[3L]] == 'key') {
                ## ASK Anqi about this
                set_slice_data(linC, linR)  ## TODO
            } else if(is.numeric(linR$data) || is.integer(linR$data))
                linC$dense_data <- format_matrix(linR$data, 'scalar')
            else if(linR$data$class == 'LinOp' && linR$data$type == 'scalar_const')
                linC$dense_data <- format_matrix(linR$data$data, 'scalar')
            else
                set_matrix_data(linC, linR)
        }
    }

    root_linC
}

