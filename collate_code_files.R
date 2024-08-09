## This script collates appropriate R files in the `rsrc_tree` directory in alphabetical order
## and puts them into one giant file in the R directory.
f <- list.files("rsrc_tree", "*.R$", full.names = TRUE, recursive = TRUE)
