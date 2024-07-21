## CVXPY SOURCE: cvxpy/atoms/affine/bmat.py

Bmat <- function(block_lists) {
  row_blocks <- lapply(block_lists, function(blocks) { do.call("HStack", blocks) })
  do.call("VStack", row_blocks)
}
