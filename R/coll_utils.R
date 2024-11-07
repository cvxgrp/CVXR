#### #' Create a vec using a list
#### #' @return a list that provides a vector as well as a deque
make_vec <- function() {
  vec <- vector(mode = "list", length = 29L) ## initial length of 29!
  n <- 0L
  push_back <- function(what) { n <<- n + 1L; vec[[n]] <<- what; invisible(what) }
  push_front <- function(what) { vec <<- append(list(what), vec); n <<- n + 1L; invisible(what) }
  pop_back <- function() { result <- vec[[n]]; vec[[n]] <<- NULL; n <<- n - 1L; result }
  pop_front <- function() { result <- vec[[1L]]; vec <<- vec[-1]; n <<- n - 1L; result }
  get_list <- function() vec[seq_len(n)]
  size <- function() n
  list(push_back = push_back, push_front = push_front, pop_back = pop_back, pop_front = pop_front,
       get_list = get_list, size = size)
}
