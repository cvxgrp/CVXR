## These tests are problematic on CRAN because they
## actually install RMOSEK and Rcplex on debian and
## they never tell us what IBM cplex is being used.
## So there is no real way to troubleshoot, short of
## installing everything ourselves.

TOL <- 1e-4

x_bool <- Variable(boolean=TRUE)
y_int <- Variable(integer=TRUE)
A_bool <- Variable(3, 2, boolean=TRUE)
B_int <- Variable(2, 3, integer=TRUE)
## Check for all installed QP solvers

MIP_SOLVERS <- c("ECOS_BB", "GUROBI", "MOSEK")

solvers <- intersect(MIP_SOLVERS, installed_solvers())

bool_prob  <- function(solver) {
    test_that("Test Boolean problems", {
        #skip_on_cran()
        ## Bool in objective
        obj <- Minimize((x_bool - 0.2)^2)
        p <- Problem(obj, list())
        result <- solve(p, solver = solver, verbose = TRUE)
        expect_equal(result$value, 0.04, tolerance = TOL)
        expect_equal(result$getValue(x_bool), 0, tolerance = TOL)

                                        # Bool in constraint
        t <- Variable()
        obj <- Minimize(t)
        p <- Problem(obj, list(x_bool^2 <= t))
        result <- solve(p, solver = solver, verbose = TRUE)
        expect_equal(result$value, 0, tolerance = TOL)
        expect_equal(result$getValue(x_bool), 0, tolerance = 1e-4)

                                        # Matrix Bool in objective
        C <- cbind(c(0,1,0), c(1,1,1))
        obj <- Minimize(sum_squares(A_bool - C))
        p <- Problem(obj, list())
        result <- solve(p, solver = solver, verbose = TRUE)
        expect_equal(result$value, 0, tolerance = TOL)
        expect_equal(result$getValue(A_bool), C, tolerance = 1e-4)

                                        # Matrix Bool in constraint
        t <- Variable()
        obj <- Minimize(t)
        p <- Problem(obj, list(sum_squares(A_bool - C) <= t))
        result <- solve(p, solver = solver, verbose = TRUE)
        expect_equal(result$value, 0, tolerance = TOL)
        expect_equal(result$getValue(A_bool), C, tolerance = 1e-4)
    })
}

int_prob  <- function(solver) {
    test_that("Test Integer problems", {
        #skip_on_cran()

        ## Int in objective
        obj <- Minimize((y_int - 0.2)^2)
        p <- Problem(obj, list())
        result <- solve(p, solver = solver, verbose = TRUE)
        expect_equal(result$value, 0.04, tolerance = TOL)
        expect_equal(result$getValue(y_int), 0, tolerance = TOL)

        ## Infeasible Int prblem
        obj <- Minimize(0)
        p <- Problem(obj, list(y_int == 0.5))
        result <- solve(p, solver = solver, verbose = TRUE)
        expect_true(result$status %in% CVXR:::INF_OR_UNB)
    })
}

int_socp  <- function(solver) {
    test_that("Test SOCP problems", {
        #skip_on_cran()

        ## Int in objective
        t  <- Variable()
        obj <- Minimize(t)
        p <- Problem(obj, list(square(y_int - 0.2) <= t))
        result <- solve(p, solver = solver, verbose = TRUE)
        expect_equal(result$value, 0.04, tolerance = TOL)
        expect_equal(result$getValue(y_int), 0, tolerance = TOL)
    })
}

bool_socp  <- function(solver) {
    test_that("Test Bool SOCP problems", {
        #skip_on_cran()

        ## Int in objective
        t  <- Variable()
        obj <- Minimize(t)
        p <- Problem(obj, list(square(x_bool - 0.2) <= t))
        result <- solve(p, solver = solver, verbose = TRUE)
        expect_equal(result$value, 0.04, tolerance = TOL)
        expect_equal(result$getValue(x_bool), 0, tolerance = TOL)
    })
}


test_all_solvers <- function() {
    for (solver in solvers) {
        bool_prob(solver)
        int_prob(solver)
        bool_socp(solver)
        int_socp(solver)
    }
}

test_all_solvers()
