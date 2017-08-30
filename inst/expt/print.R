
linop2str = function(linop) {
    class <- linop$class
    type <- linop$type
    size <- paste(linop$size, collapse=", ")
    ## args <- if (is.null(linop$args)) "NULL" else linop2str(linop$args)
    data <- linop$data
    ## list(class, type, size, args, data)
    cat(class, type, size, data, "\n")
    sprintf("%s(type=%s, size=[%s], data=[%s])",
            linop$class,
            linop$type,
            paste(linop$size, collapse=", "),
            ##if (is.null(linop$args)) "NULL" else linopList2str(linop$args),
            if (is.null(linop$data) || length(linop$data == 0)) "" else linop$data)
}


linopList2str = function(linopList) {
    result <- sapply(linopList, linop2str)
    result <- paste(result, collapse = ", ")
    sprintf("[ %s ]", result)
}

