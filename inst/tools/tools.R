## Crude tools for my use; will polish up later (BN)
makeR6Class <- function(className, publicFields=NULL, publicMethods=NULL, filename=paste0(className, "-R6.R")) {
    template <- readLines("template.R")
    template <- gsub("<CLASSNAME>", className, template)
    activeStart <- grep("<ACTIVE_START>", template)
    activeEnd <- grep("<ACTIVE_END>", template)
    publicStart <- grep("<PUBLIC_START>", template)
    publicEnd <- grep("<PUBLIC_END>", template)

    activeChunk <- template[(activeStart + 1):(activeEnd - 1)]
    publicChunk <- template[(publicStart + 1):(publicEnd - 1)]
    activeLines <- lapply(publicFields, function(x) gsub("<VAR>", x, activeChunk))
    publicLines <- lapply(publicMethods, function(x) gsub("<NAME>", x, publicChunk))
    result <- template[1:(activeStart - 1)]
    for (i in seq(activeLines)) {
        if (i > 1) {
            result <- c(result, ",")
        }
        result <- c(result, activeLines[[i]])
    }
    result <- c(result, template[(activeEnd + 1):(publicStart - 1)])
    for (i in seq(publicLines)) {
        result <- c(result, ",", publicLines[[i]])
    }
    result <- c(result, template[(publicEnd + 1):length(template)])
    writeLines(text = result, con = paste0(className, "-R6.R"))
}


makeR6Class(className = "CVXCanon.ProblemData",
            publicFields = c("V", "I", "J", "id_to_col", "const_to_row",
                             "num_constraints", "vals", "row_idxs", "col_ptrs"),
            publicMethods = c("getV", "getI", "getJ", "getConstVec"))


makeR6Class(className = "CVXCanon.Solution",
            publicFields = c("status", "optimal_value", "primal_values", "dual_values"),
            publicMethods = c())


makeR6Class(className = "CVXCanon.LinOp",
            publicFields = c("sparse_data", "dense_data", "type", "slice"),
            publicMethods = c("args_push_back", "size_push_back"))

makeR6Class(className = "CVXCanon.LinOpVector",
            publicFields = c(),
            publicMethods = c("push_back"))
