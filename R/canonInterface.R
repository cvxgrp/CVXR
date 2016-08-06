## Added format_matrix and set_matrix_data.
get_problem_matrix <- function(constrs, id_to_col = NA, constr_offsets = NA) {
    linOps <- lapply(constrs, function(constr) { constr$expr })
    lin_vec <- CVXcanon.LinOpVector$new()

    id_to_col_C <- CVXcanon.IntIntMap$new()
    if (is.na(id_to_col))
        id_to_col <- list()

    ## Loading the variable offsets from our R list into a C++ map
    for (id in names(id_to_col)) {
        col <- id_to_col[[id]]
        id_to_col_C[as.integer(id)] <- as.integer(col)
    }

    ## This array keeps variables data in scope after build_lin_op_tree returns
    tmp <- R6List$new()
    for (lin in linOps) {
        tree <- build_lin_op_tree(lin, tmp)
        tmp$append(tree)
        lin_vec.push_back(tree)
    }

    if (is.na(constr_offsets))
        problemData <- CVXcanon.build_matrix(lin_vec, id_to_col_C)
    else {
        ## Load constraint offsets into a C++ vector
        constr_offsets_C <- CVXcanon.IntVector()
        for (offset in constr_offsets)
            constr_offsets_C.push_back(as.integer(offset))
        problemData <- CVXcanon.build_matrix(lin_vec, id_to_col_C, constr_offsets_C)
    }

    ## Unpacking
    V <- problemData.getV(length(problemData.V))
    I <- problemData.getI(length(problemData.I))
    J <- problemData.getJ(length(problemData.J))
    const_vec <- problemData.getConstVec(length(problemData.const_vec))

    list(V = V, I = I, J = J, const_vec = const_vec.reshape(-1, 1))
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

## What does this do?
set_slice_data <- function(linC, linR) {  ## What does this do?
    for (i in seq_int(length(linR$data))) {
        sl <- linR$data[[i]]
        #$vec <- CVXcanon.IntVector()
        vec <- list()

        ## In R this is just a list of int vectors.

        start_idx <- 0
        if (!is.na(sl$start_idx))
            start_idx <- sl$start_idx

        stop_idx <- linR$args[[1]]$size[i]
        if (!is.na(sl$stop_idx))
            stop_idx <- sl$stop_idx

        step <- 1
        if(!is.na(sl$step))
            step <- sl$step

        ## handle [::-1] case
        if(step < 0 && is.na(sl$start_idx) && is.na(sl$stop_idx)) {
            tmp <- start
            start_idx <- stop_idx - 1
            stop_idx <- tmp
        }

        for(var in c(start_idx, stop_idx, step))
            vec.push_back(var)

        linC.slice.push_back(vec)
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
            if(length(linR$data) == 2 && is(linR$data[1], 'slice'))
                ## ASK Anqi about this
                set_slice_data(linC, linR)  ## TODO
            else if(is.numeric(linR$data) || is.integer(linR$data))
                linC$dense_data <- format_matrix(linR$data, 'scalar')
            else if(linR$data$class == 'LinOp' && linR$data$type == 'scalar_const')
                linC$dense_data <- format_matrix(linR$data$data, 'scalar')
            else
                set_matrix_data(linC, linR)
        }
    }

    root_linC
}

get_constraint_node <- function(c, tmp) {
    ## c is the constraint, tmp is what is returned
    root <- CVXcanon.LinOp$new() ## create a C linop
    if(is.list(c)) {
        c_size <- c$size
        c_constr_id <- c$constr_id
    } else {
        c_size <- size(c)
        c_constr_id <- c@constr_id
    }
    root$size_push_back(c_size[1])
    root$size_push_back(c_size[2])

    ## add ID as dense_data
    root$dense_data <- format_matrix(c_constr_id, "scalar")

    if(is.list(c) && c$class == "LinEqConstr") {
        root$type <- LINOP_TYPES["EQ"]
        expr <- build_lin_op_tree(c$expr, tmp)
        tmp$append(expr)
        root$args_push_back(expr)
    } else if(is.list(c) && c$class == "LinLeqConstr") {
        root$type <- LINOP_TYPES["LEQ"]
        expr <- build_lin_op_tree(c$expr, tmp)
        tmp$append(expr)
        root$args_push_back(expr)
    } else if(is(c, "SOC")) {
        root$type <- LINOP_TYPES["SOC"]
        tt <- build_lin_op_tree(c@t, tmp)
        tmp$append(tt)
        root$args_push_back(tt)
        for(elem in c@x_elems) {
            x_elem <- build_lin_op_tree(elem, tmp)
            tmp$append(x_elem)
            root$args_push_back(x_elem)
        }
    } else if(is(c, "ExpCone")) {
        root$type <- LINOP_TYPES["EXP"]
        x <- build_lin_op_tree(c@x, tmp)
        y <- build_lin_op_tree(c@y, tmp)
        z <- build_lin_op_tree(c@z, tmp)
        root$args_push_back(x)
        root$args_push_back(y)
        root$args_push_back(z)
        tmp$append(list(x, y, z))
    } else if(is(c, "SDP"))
        stop("Unimplemented: SDP")
    else
        stop("Undefined constraint type ", class(c))
    return(root)
}

solve <- function(sense, objective, constraints, verbose, solver_options = ecos.control()) {
    ## This array keeps variables data in scope after build_lin_op_tree returns
    tmp <- R6List$new()

    C_objective <- build_lin_op_tree(objective, tmp)

    C_constraints <- CVXcanon.LinOpVector$new()
    for(constr in constraints) {
        ## Gets the constraint node
        root <- get_constraint_node(constr, tmp)
        tmp$append(root)
        C_constraints$push_back(root)
    }

    cat("Objective:\n")
    print(C_objective$toString())

    cat("Constraints:\n")
    print(C_constraints$toString())

    if(sense == "Minimize")
        C_sense <- 0L ## CVXcanon.MINIMIZE
    else if(sense == "Maximize")
        C_sense <- 1L ##CVXcanon.MAXIMIZE
    else
        stop("Unimplemented")

    C_opts <- solver_options
    C_opts['verbose'] <- 1L

    C_opts <- list(verbose = 1L)
    cvxCanon <- CVXcanon$new()
    solution <- cvxCanon$solve(C_sense, C_objective, C_constraints, C_opts)

    ## print(paste("CVXcanon optimal value", solution$optimal_value))

    return(solution)
}
