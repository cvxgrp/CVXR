## How to use `pickle` for R/Python

```{r}
#' Pickle an R object to load in python
#' @param x an R object that is transferrable to python, e.g. list of matrices, vectors, etc.
#' @param fname the name to use for the pickle file
#' @examples
#' pickle(list(a = 1.0, b = 2L, x = matrix(1:20, nrow = 5), y = rnorm(5)), "foo.pickle")
pickle <- function(x, fname) {
  io <- reticulate::import("io")
  f <- io$open(fname, mode = "wb")
  pk <- reticulate::import("pickle")
  pk$dump(x, f)
  f$close()
  ## Then, on the python side
  ## import pickle as pk
  ## with open("foo.pickle", "rb") as f:
  ##    l = pk.load(f)
  ## Now l is a dict on the python side with components l['a'], l['x'] etc.
}
#' Unpickle a Python dict to load in R
#' @param fpath the file path for the pickle file containing a python dict of objects transferrable to R, such as matrices, vectors, etc.
#' @examples
#' l <- unpickle("./foo.pickle")  ## l is now a named list
#' 
unpickle <- function(fname) {
  io <- reticulate::import("io")
  py_code <- c(
    "import pickle",
    "with open('%s', 'rb') as file:",
    "    result = pickle.load(file)"
  )
  py_code[2] <- sprintf(py_code[2], fname)
  tmp_file <- tempfile(fileext = ".py")
  writeLines(py_code, tmp_file)
  py_run_file(tmp_file)
  py$result
}
```
