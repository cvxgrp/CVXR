LINOP_TYPES <- c(VARIABLE = "VARIABLE",
                 PROMOTE = "PROMOTE",
                 MUL_ELEM = "__ABCDE__",
                 RMUL = "__FGHIJ__",
                 MUL = "MUL",
                 DIV = "DIV",
                 SUM_ENTRIES = "__KLMNO__",
                 SUM = "SUM",
                 NEG = "NEG",
                 INDEX = "INDEX",
                 TRANSPOSE = "TRANSPOSE",
                 TRACE = "TRACE",
                 RESHAPE = "RESHAPE",
                 DIAG_VEC = "DIAG_VEC",
                 DIAG_MAT = "DIAG_MAT",
                 UPPER_TRI = "UPPER_TRI",
                 CONV = "CONV",
                 HSTACK = "HSTACK",
                 VSTACK = "VSTACK",
                 SCALAR_CONST = "SCALAR_CONST",
                 DENSE_CONST = "DENSE_CONST",
                 SPARSE_CONST = "SPARSE_CONST",
                 NO_OP = "NO_OP",
                 KRON = "KRON",
                 LEQ = "__PQRST__",
                 EQ = "EQ",
                 SOC = "SOC",
                 EXP = "EXP",
                 SDP = "SDP")

files <- list.files("../../R", ".R$", full.names=TRUE)

dir <- "tmp"

for (x in files) {
    print(x)
    lines <- readLines(x)
    for (y in names(LINOP_TYPES)) {
        lines <- gsub(y, paste('LINOP_TYPES["', LINOP_TYPES[y], '"]', sep=""), lines)
    }
    lines <- gsub("__ABCDE__", "MUL_ELEM", lines)
    lines <- gsub("__FGHIJ__", "RMUL", lines)
    lines <- gsub("__KLMNO__", "SUM_ENTRIES", lines)
    lines <- gsub("__PQRST__", "LEQ", lines)

    writeLines(lines, con= paste(dir, basename(x), sep="/"))
}
