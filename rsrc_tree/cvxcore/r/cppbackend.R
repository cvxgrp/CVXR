## CVXPY SOURCE: cvxpy/cvxcore/python/cppbackend.py

build_matrix <- function(id_to_col,
                         param_to_size,
                         param_to_col,
                         var_length,
                         constr_length,
                         linOps) {
  # return a CSC sparse matrix
  lin_vec <- cvxcore.ConstLinOpVector$new()
  id_to_col_C <- id_to_col # No need for conversion since it is already a named vector, with ints for names
  param_to_size_C <- param_to_size # No need for conversion since it is already a named vector, with ints for names

  # dict to memoize construction of C++ linOps, and to keep Python references
  # to them to prevent their deletion
  linPy_to_linC <- list()
  for (lin in linOps) {
    build_lin_op_tree(lin, linPy_to_linC)
    tree <- linPy_to_linC[lin]
    lin_vec$push_back(tree)
  }

    problemData = cvxcore.build_matrix(
        lin_vec, int(var_length), id_to_col_C, param_to_size_C, s.get_num_threads()
    )

    # Populate tensors with info from problemData.
    tensor_V = {}
    tensor_I = {}
    tensor_J = {}
    for param_id, size in param_to_size.items():
        tensor_V[param_id] = []
        tensor_I[param_id] = []
        tensor_J[param_id] = []
        problemData.param_id = param_id
        for i in range(size):
            problemData.vec_idx = i
            prob_len = problemData.getLen()
            tensor_V[param_id].append(problemData.getV(prob_len))
            tensor_I[param_id].append(problemData.getI(prob_len))
            tensor_J[param_id].append(problemData.getJ(prob_len))

    # Reduce tensors to a single sparse CSR matrix.
    V = []
    I = []
    J = []
    # one of the 'parameters' in param_to_col is a constant scalar offset,
    # hence 'plus_one'
    param_size_plus_one = 0
    for param_id, col in param_to_col.items():
        size = param_to_size[param_id]
        param_size_plus_one += size
        for i in range(size):
            V.append(tensor_V[param_id][i])
            I.append(tensor_I[param_id][i] + tensor_J[param_id][i] * constr_length)
            J.append(tensor_J[param_id][i] * 0 + (i + col))
    V = np.concatenate(V)
    I = np.concatenate(I)
    J = np.concatenate(J)

    output_shape = (
        np.int64(constr_length) * np.int64(var_length + 1),
        param_size_plus_one,
    )
    A = sp.csc_matrix((V, (I, J)), shape=output_shape)
    return A

}


#' Set numerical data fields in linC
#'
#' This function sets numerical data fields in the linC object based on the
#' linR object's data.
#'
#' @param linC An object in which the data fields will be set.
#' @param linR An object that contains the data to be set in linC.
#' @return None
#' @export
set_linC_data <- function(linC, linR) {
  ## Assert that linR$data is not NULL
  if (is.null(linR$data)) stop(sprintf("linR(%s) data is NULL!", linR$uuid))
  
  ## If linR$data is a tuple and the first element is a slice, call set_slice_data
  ## We have to detect a "slice" as a sequence of three items
  ## RIGHT NOW I AM JUST USING WHAT WE DID BEFORE, NOT SURE IT IS CORRECT!
  ## ASK ANQI
  if (length(linR$data) == 3L && linR$data[[3L]] == 'key') {
    ## ASK Anqi about this
    set_slice_data(linC, linR)  ## TODO
    ## If linR$data is a float or an integer, set dense data and data dimension
  } else if (is.numeric(linR$data) || is.integer(linR$data)) {
    linC$set_dense_data(format_matrix(linR$data, format='scalar'))
    linC$set_data_ndim(0)
  # Otherwise, call set_matrix_data
  } else {
    set_matrix_data(linC, linR)
  }
}

#' Construct a C++ LinOp corresponding to LinPy.
#'
#' Children of linR are retrieved from linR_to_linC.
#'
#' @param linR A LinPy object.
#' @param linR_to_linC A mapping from LinPy objects to LinC objects.
make_linC_from_linR <- function(linR, linR_to_linC) {
  if (! (linR$uuid %in% names(linR_to_linC))) {
    typ <- linR$type
    dim <- linR$dim
    lin_args_vec <- make_ConstLinOpVector()
    for (arg in linR$args) {
      lin_args_vec$push_back(linR_to_linC[[arg]])
    }
    
    linC <- make_cvxcore_LinOp(type = typ, dim = dim, args = lin_args_vec)
    linR_to_linC[[linR$uuid]] <- linC
    
    if (!is.null(linR$data)) {
      if (inherits(linR$data, "LinOp")) {
        linR_data <- linR$data
        linC_data <- linR_to_linC[[linR_data$uuid]]
        linC$set_linOp_data(linC_data)
        linC$set_data_ndim(length(linR_data$dim))
      } else {
        set_linC_data(linC, linR)
      }
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

        ## stop_idx <- linR$args[[1]]$dim[i] - 1L
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
        ## if (length(sl) == 1L) {
        ##     vec <- c(sl - 1L, sl, 1L)
        ## } else if (length(sl) == 2L) {
        ##     vec <- c(sl[1L] - 1L, sl[2L], 1L)  # Using zero-based indexing, and step assumed to be 1.
        ## } else {
        ##     r <- range(sl)
        ##     vec <- c(r[1L] - 1L, r[2L], 1L)
        ## }

        ##vec <- c(sl, 1L)  # Using 1-based indexing, and step assumed to be 1.
        linC$slice_push_back(sl - 1L) ## Make indices zero-based for C++
    }
}

#' Construct C++ LinOp tree from Python LinOp tree.
#'
#' Constructed C++ linOps are stored in the linR_to_linC list,
#' which maps Python linOps to their corresponding C++ linOps.
#'
#' @param root_linR The root of the Python LinOp tree.
#' @param linR_to_linC A list for memoizing construction and storing
#'   the C++ LinOps.
build_lin_op_tree <- function(root_linR, linR_to_linC) {
  bfs_stack <- list(root_linR)
  post_order_stack <- list()
  
  while (length(bfs_stack) > 0) {
    n <- length(bfs_stack)
    linR <- bfs_stack[[n]]
    bfs_stack[[n]] <- NULL ## pop off the last element
    if (linR$uuid %in% names(linR_to_linC)) {
      n <- length(pos_order_stack)
      post_order_stack[[n + 1L]] <- linR
      for (arg in linR$args) {
        n <- length(bfs_stack)
        bfs_stack[[n + 1L]] <- arg
      }
      if (inherits(linR$data, "LinOp")) {
        n <- length(bfs_stack)
        bfs_stack[[n + 1L]] <- linR$data
      }
    }
  }
  
  while (length(post_order_stack) > 0) {
    n <- length(post_order_stack)
    linR <- post_order_stack[[n]]
    post_order_stack[[n]] <- NULL ## delete linR
    make_linC_from_linR(linR, linR_to_linC)
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

format_matrix <- function(matrix, format='dense') {
    ## Returns the matrix in the appropriate form,
    ## so that it can be efficiently loaded with our swig wrapper

    ## TODO: Should we convert bigq/bigz values? What if it's a sparse matrix?
    if(is.bigq(matrix) || is.bigz(matrix)) {
        matdbl <- matrix(sapply(matrix, as.double))
        dim(matdbl) <- dim(matrix)
        matrix <- matdbl
    }

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




