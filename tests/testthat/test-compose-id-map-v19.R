## @cvxpy tests/unit/reductions/test_chain.py
## Chain.compose_var_id_map / compose_param_id_map (#3147 part A surface).
##
## CVXR keeps the shared-id design (part B deferred, ADR D_19.5), so on a real
## chain these compose to list(). But the composition MACHINERY is a faithful
## port of chain.py:_compose_id_map -- these tests prove it (a) composes 1->many
## correctly on synthetic maps, and (b) the consumer idiom falls back to the
## original id under shared-id. This is the "cheap insurance" that lets a future
## CVXPY consumer port work verbatim without flipping the producers.

ns <- asNamespace("CVXR")
.compose_id_map <- get(".compose_id_map", ns)

test_that(".compose_id_map composes a chain of 1->many maps (chain.py example)", {
  ## reduction i: A -> [A']; reduction j>i: A' -> [A'', A''']  =>  A -> [A'', A''']
  step1 <- list("A" = "A'")
  step2 <- list("A'" = c("A''", "A'''"))
  out <- .compose_id_map(list(step1, step2))
  expect_identical(out[["A"]], c("A''", "A'''"))
})

test_that(".compose_id_map carries unmapped ids through and adds new ones", {
  ## B is untouched by step1, introduced fresh-mapped by step2; C only in step1.
  step1 <- list("C" = "C'")
  step2 <- list("B" = "B'")
  out <- .compose_id_map(list(step1, step2))
  expect_identical(out[["C"]], "C'")   # forwarded (step2 doesn't touch C')
  expect_identical(out[["B"]], "B'")   # newly introduced by step2
})

test_that(".compose_id_map over empty step maps is empty (shared-id case)", {
  expect_identical(.compose_id_map(list()), list())
  expect_identical(.compose_id_map(list(list(), list())), list())
})

test_that("compose_var_id_map/compose_param_id_map are list() on a shared-id chain", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - 1)), list(x >= 0))
  chain_obj <- get("construct_solving_chain", ns)(prob, solver = "CLARABEL")
  cvar <- get("compose_var_id_map", ns)(chain_obj)
  cpar <- get("compose_param_id_map", ns)(chain_obj)
  expect_identical(cvar, list())
  expect_identical(cpar, list())

  ## consumer idiom: compose_var_id_map(chain)[[id]] %||% id  -> original id
  vid <- as.character(x@id)
  resolved <- cvar[[vid]]
  if (is.null(resolved)) resolved <- vid
  expect_identical(resolved, vid)
})
