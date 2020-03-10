## These tests were moved here because scs causes errors on R-devel on travis
## Also because they go over the time limit
context("test-examples")
TOL <- 1e-6

test_that("Test Chebyshev center", {
    skip_on_cran()
    # Find the largest Euclidean ball in the polyhedron.
    # Goal is to find the largest Euclidean ball (i.e. its center and radius)
    # that lies in a polyhedron described by linear inequalities in this fashion:
    # P = {x : a_i'*x <= b_i, i=1,...,m} where x is in R^2.

    # Generate the input data.
    a1 <- matrix(c(2,1))
    a2 <- matrix(c(2,-1))
    a3 <- matrix(c(-1,2))
    a4 <- matrix(c(-1,-2))
    b <- matrix(1, nrow = 4)

    # Create and solve the model.
    r <- Variable(name = "r")
    x_c <- Variable(2, name = "x_c")
    obj <- Maximize(r)
    constraints <- list(   # TODO: Have atoms compute values for constants.
        t(a1) %*% x_c + norm(a1, "F")*r <= b[1],
        t(a2) %*% x_c + norm(a2, "F")*r <= b[2],
        t(a3) %*% x_c + norm(a3, "F")*r <= b[3],
        t(a4) %*% x_c + norm(a4, "F")*r <= b[4]
    )

    p <- Problem(obj, constraints)
    result <- solve(p)
    expect_equal(result$value, 0.447214, tolerance = TOL)
    expect_equal(result$getValue(r), result$value, tolerance = TOL)
    expect_equal(result$getValue(x_c), matrix(c(0,0)), tolerance = TOL)
})

test_that("Test issue with scalars", {
    skip_on_cran()

    n <- 6
    eps <- 1e-6
    set.seed(10)
    P0 <- matrix(stats::rnorm(n*n), ncol = n)
    eye <- diag(n)
    P0 <- t(P0) %*% P0 + eps * eye

    print(P0)

    P1 <- matrix(stats::rnorm(n*n), nrow = n, ncol = n)
    P1 <- t(P1) %*% P1
    P2 <- matrix(stats::rnorm(n*n), nrow = n, ncol = n)
    P2 <- t(P2) %*% P2
    P3 <- matrix(stats::rnorm(n*n), nrow = n, ncol = n)
    P3 <- t(P3) %*% P3

    q0 <- matrix(stats::rnorm(n), ncol = 1)
    q1 <- matrix(stats::rnorm(n), ncol = 1)
    q2 <- matrix(stats::rnorm(n), ncol = 1)
    q3 <- matrix(stats::rnorm(n), ncol = 1)

    r0 <- matrix(stats::rnorm(1), ncol = 1)
    r1 <- matrix(stats::rnorm(1), ncol = 1)
    r2 <- matrix(stats::rnorm(1), ncol = 1)
    r3 <- matrix(stats::rnorm(1), ncol = 1)

    slack <- Variable()

    # Form the problem.
    x <- Variable(n)
    objective <- Minimize(0.5*quad_form(x, P0) + t(q0) %*% x + r0 + slack)
    constraints <- list(0.5*quad_form(x, P1) + t(q1) %*% x + r1 <= slack,
                        0.5*quad_form(x, P2) + t(q2) %*% x + r2 <= slack,
                        0.5*quad_form(x, P3) + t(q3) %*% x + r3 <= slack
    )

    # We now find the primal result and compare it to the dual result to check
    # if strong duality holds, i.e., the duality gap is effectively zero.
    p <- Problem(objective, constraints)
    result <- solve(p)

    # Note that since our data is random, we may need to run this program multiple times
    # to get a feasible primal. When feasible, we can print out the following values:
    print(result$getValue(x))   # solution
    lam1 <-  result$getDualValue(constraints[[1]])
    lam2 <-  result$getDualValue(constraints[[2]])
    lam3 <-  result$getDualValue(constraints[[3]])
    print(class(lam1))

    P_lam <- P0 + lam1*P1 + lam2*P2 + lam3*P3
    q_lam <- q0 + lam1*q1 + lam2*q2 + lam3*q3
    r_lam <- r0 + lam1*r1 + lam2*r2 + lam3*r3
    dual_result <- -0.5*t(q_lam) %*% P_lam %*% q_lam + r_lam
    print(dim(dual_result))
    expect_equal(CVXR:::intf_dim(dual_result), c(1,1))
})

test_that("Test examples from the README", {
    skip_on_cran()

    ## Problem data
    m <- 30
    n <- 20
    A <- matrix(stats::rnorm(m*n), nrow = m, ncol = n)
    b <- matrix(stats::rnorm(m), nrow = m, ncol = 1)

    ## Construct the problem
    x <- Variable(n)
    objective <- Minimize(sum((A %*% x - b)^2))
    constraints <- list(x >= 0, x <= 1)
    p <- Problem(objective, constraints)

    ## The optimal objective is returned by solve(p)
    result <- solve(p)
    ## The optimal value for x is stored in result$getValue(x)
    print(result$getValue(x))
    ## The optimal Lagrange multiplier for a constraint is stored in result$getDualValue(constraints[[1]])
    print(result$getDualValue(constraints[[1]]))

    ###########################################
    ## Scalar variable
    a <- Variable()

    ## Column vector variable of length 5
    x <- Variable(5)

    ## Matrix variable with 4 rows and 7 columns
    A <- Variable(4, 7)

    ###########################################
    ## Positive scalar parameter
    m <- Parameter(nonneg = TRUE)

    ## Column vector parameter with unknown sign (by default)
    c <- Parameter(5)

    ## Matrix parameter with negative entries
    G <- Parameter(4, 7, nonpos = TRUE)

    ## Assigns a constant value to G
    value(G) <- -matrix(1, nrow = 4, ncol = 7)

    # Raises an error for assigning a value with invalid sign
    expect_error(value(G) <- matrix(1, nrow = 4, ncol = 7), "Value must be nonpositive")

    ###########################################
    a <- Variable()
    x <- Variable(5)

    ## expr is an Expression object after each assignment
    expr <- 2*x
    expr <- expr - a
    expr <- sum(expr) + norm2(x)

    ###########################################
    ## Problem data
    n <- 10
    m <- 5
    A <- matrix(stats::rnorm(n*m), nrow = n, ncol = m)
    b <- matrix(stats::rnorm(n), nrow = n, ncol = 1)
    gamma <- Parameter(nonneg = TRUE)

    ## Construct the problem
    x <- Variable(m)
    loss <- sum((A %*% x - b)^2)

    get_x <- function(gamma_value) {
        value(gamma) <- gamma_value
        objective <- Minimize(loss + gamma*norm1(x))
        # objective <- Minimize(loss + gamma_value*norm1(x))
        p <- Problem(objective)

        result <- solve(p)
        result$getValue(x)
    }

    gammas <- 10^seq(-1, 2, length.out = 2)

    ## Serial computation
    x_values <- sapply(gammas, get_x)

    ###########################################
    n <- 10

    mu <- matrix(stats::rnorm(n), nrow = 1, ncol = n)
    sigma <- matrix(stats::rnorm(n*n), nrow = n, ncol = n)
    sigma <- t(sigma) %*% sigma
    gamma <- Parameter(nonneg = TRUE)
    value(gamma) <- 1
    x <- Variable(n)

    ## Constants:
    ## mu is the vector of expected returns.
    ## sigma is the covariance matrix.
    ## gamma is a Parameter that trades off risk and return.

    ## Variables:
    ## x is a vector of stock holdings as fractions of total assets.

    expected_return <- mu %*% x
    risk <- quad_form(x, sigma)

    objective <- Maximize(expected_return - gamma*risk)
    p <- Problem(objective, list(sum(x) == 1))
    result <- solve(p)

    ## The optimal expected return
    print(result$getValue(expected_return))

    ## The optimal risk
    print(result$getValue(risk))

    ###########################################
    N <- 50
    M <- 40
    n <- 10
    data1 <- lapply(1:N, function(i) { list(1, matrix(stats::rnorm(n, mean = 1.0, sd = 2.0))) })
    data2 <- lapply(1:M, function(i) { list(-1, matrix(stats::rnorm(n, mean = -1.0, sd = 2.0))) })
    data <- c(data1, data2)

    ## Construct problem
    gamma <- Parameter(nonneg = TRUE)
    value(gamma) <- 0.1
    ## 'a' is a variable constrained to have at most 6 nonzero entries
    a <- Variable(n)
    b <- Variable()

    slack <- lapply(data, function(x) {
        label <- x[[1]]
        sample <- x[[2]]
        pos(1 - label*(t(sample) %*% a - b))
    })
    objective <- Minimize(norm2(a) + gamma*Reduce("+", slack))
    p <- Problem(objective)
    result <- solve(p)

    ## Count misclassifications
    errors <- 0
    for(v in data) {
        label <- v[[1]]
        sample <- v[[2]]
        if(label * result$getValue(t(sample) %*% a - b) < 0)
            errors <- errors + 1
    }

    print(paste(errors, "misclassifications"))
    print(result$getValue(a))
    print(result$getValue(b))
})

test_that("Test advanced tutorial 1", {
    skip_on_cran()

    # Solving a problem with different solvers.
    x <- Variable(2)
    obj <- Minimize(x[1] + norm1(x))
    constraints <- list(x >= 2)
    prob <- Problem(obj, constraints)

    # Solve with ECOS.
    result <- solve(prob, solver = "ECOS")
    print(paste("Optimal value with ECOS:", result$value))
    expect_equal(result$value, 6, tolerance = TOL)

    # Solve with ECOS_BB.
    result <- solve(prob, solver = "ECOS_BB")
    print(paste("Optimal value with ECOS_BB:", result$value))
    expect_equal(result$value, 6, tolerance = TOL)

    # Solve with CVXOPT.
    if("CVXOPT" %in% installed_solvers()) {
        result <- solve(prob, solver = "CVXOPT")
        print(paste("Optimal value with CVXOPT:", result$value))
        expect_equal(result$value, 6, tolerance = TOL)
    }

    # Solve with SCS.
    result <- solve(prob, solver = "SCS")
    print(paste("Optimal value with SCS:", result$value))
    expect_equal(result$value, 6, tolerance = 1e-2)

    # Solve with CPLEX.
    if("CPLEX" %in% installed_solvers()) {
        result <- solve(prob, solver = "CPLEX")
        print(paste("Optimal value with CPLEX:", result$value))
        expect_equal(result$value, 6, tolerance = TOL)
    }

    if("GLPK" %in% installed_solvers()) {
        # Solve with GLPK.
        result <- solve(prob, solver = "GLPK")
        print(paste("Optimal value with GLPK:", result$value))
        expect_equal(result$value, 6, tolerance = TOL)

        # Solve with GLPK_MI.
        result <- solve(prob, solver = "GLPK_MI")
        print(paste("Optimal value with GLPK_MI:", result$value))
        expect_equal(result$value, 6, tolerance = TOL)
    }

    # Solve with MOSEK.
    if("MOSEK" %in% installed_solvers()) {
        result <- solve(prob, solver = "MOSEK")
        print(paste("Optimal value with MOSEK:", result$value))
        expect_equal(result$value, 6, tolerance = TOL)
    }

    # Solve with GUROBI.
    if("GUROBI" %in% installed_solvers()) {
        result <- solve(prob, solver = "GUROBI")
        print(paste("Optimal value with GUROBI:", result$value))
        expect_equal(result$value, 6, tolerance = TOL)
    }

    print(installed_solvers())
})

test_that("Test log-determinant", {
    skip_on_cran()

    # Generate data
    x <- cbind(c(0.55, 0.0),
               c(0.25, 0.35),
               c(-0.2, 0.2),
               c(-0.25, -0.1),
               c(0.0, -0.3),
               c(0.4, -0.2))
    n <- nrow(x)
    m <- ncol(x)

    # Create and solve the model
    A <- Variable(n,n)
    b <- Variable(n)
    obj <- Maximize(log_det(A))
    constraints <- lapply(seq_len(m), function(i) { norm2(A %*% x[,i] + b) <= 1 })
    p <- Problem(obj, constraints)
    result <- solve(p)
    expect_equal(result$value, 1.9746, tolerance = 1e-2)
})

test_that("Test portfolio problem that caused DCPAttr errors", {
    skip_on_cran()

    library(Matrix)
    set.seed(5)
    n <- 100   # 10000
    m <- 10   # 100

    Fmat <- rsparsematrix(m, n, density = 0.01)
    Fmat@x <- rep(1, length(Fmat@x))
    D <- sparseMatrix(i = seq_len(n), j = seq_len(n), x = rep(1,n))
    D@x <- rnorm(length(D@x))^2
    Z <- matrix(rnorm(m), nrow = m, ncol = 1)
    Z <- Z %*% t(Z)

    x <- Variable(n)
    y <- Fmat %*% x

    # DCPAttr causes error because not all the curvature matrices are
    # reduced to constants when an atom is scalar.
    norm2(D %*% x)^2 + (Z %*% y)^2
})

test_that("Test examples from CVXR introduction", {
    skip_on_cran()

    m <- 30
    n <- 20
    set.seed(1)
    A <- matrix(rnorm(m*n), nrow = m, ncol = n)
    b <- matrix(rnorm(m), nrow = m, ncol = 1)

    # Construct the problem.
    x <- Variable(n)
    objective <- Minimize(sum((A %*% x - b)^2))
    constraints <- list(x >= 0, x <= 1)
    prob <- Problem(objective, constraints)

    # The optimal objective is returned by solve(prob).
    result <- solve(prob)
    # The optimal value for x is stored in result$getValue(x)
    print(result$getValue(x))
    # The optimal Lagrange multiplier for a constraint is stored in result$getDualValue(constraint[[1]])
    print(result$getDualValue(prob@constraints[[1]]))

    ########################################

    # Create the two scalar variables.
    x <- Variable()
    y <- Variable()

    # Create two constraints.
    constraints <- list(x + y == 1, x - y >= 1)

    # Form objective.
    obj <- Minimize((x - y)^2)

    # Form and solve problem.
    prob <- Problem(obj, constraints)
    result <- solve(prob)   # Returns the optimal value.
    print(paste("status:", result$status))
    print(paste("optimal value:", result$value))
    print(paste("optimal var (x):", result$getValue(x)))
    print(paste("optimal var (y):", result$getValue(y)))

    ########################################

    # Create the two scalar variables.
    x <- Variable()
    y <- Variable()

    # Create two constraints.
    constraints <- list(x + y == 1, x - y >= 1)

    # Form objective.
    obj <- Minimize((x - y)^2)

    # Form and solve problem.
    prob <- Problem(obj, constraints)
    result <- solve(prob)   # Returns the optimal value.
    print(paste("status:", result$status))
    print(paste("optimal value:", result$value))
    print(paste("optimal var (x):", result$getValue(x)))
    print(paste("optimal var (y):", result$getValue(y)))

    expect_equal(result$status, "optimal")
    expect_equal(result$value, 1.0, tolerance = TOL)
    expect_equal(result$getValue(x), 1.0, tolerance = TOL)
    expect_equal(result$getValue(y), 0, tolerance = TOL)

    ########################################

    # Replace the objective.
    prob <- Problem(Maximize(x + y), prob@constraints)
    result <- solve(prob)
    print(paste("optimal value:", result$value))

    expect_equal(result$value, 1.0, tolerance = 1e-3)

    # Replace the constraint (x + y == 1).
    constraints <- prob@constraints
    constraints[[1]] <- (x + y <= 3)
    prob <- Problem(prob@objective, constraints)
    result <- solve(prob)
    print(paste("optimal value:", result$value))

    expect_equal(result$value, 3.0, tolerance = 1e-2)

    ########################################

    x <- Variable()

    # An infeasible problem.
    prob <- Problem(Minimize(x), list(x >= 1, x <= 0))
    result <- solve(prob)
    print(paste("status:", result$status))
    print(paste("optimal value:", result$value))

    expect_equal(result$status, "infeasible")
    expect_equal(result$value, Inf)

    # An unbounded problem.
    prob <- Problem(Minimize(x))
    result <- solve(prob)
    print(paste("status:", result$status))
    print(paste("optimal value:", result$value))

    expect_equal(result$status, "unbounded")
    expect_equal(result$value, -Inf)

    ########################################

    # A scalar variable.
    Variable()

    # Column vector variable of length 5.
    x <- Variable(5)

    # Matrix variable with 4 rows and 7 columns.
    A <- Variable(4,7)

    ########################################

    m <- 10
    n <- 5
    set.seed(1)
    A <- matrix(rnorm(m*n), nrow = m, ncol = n)
    b <- matrix(rnorm(m))

    # Construct the problem.
    x <- Variable(n)
    objective <- Minimize(sum((A %*% x - b)^2))
    constraints <- list(x >= 0, x <= 1)
    prob <- Problem(objective, constraints)
    result <- solve(prob)

    print(paste("Optimal value:", result$value))
    print("Optimal var:")
    print(result$getValue(x))

    expect_equal(result$value, 7.24427715882, tolerance = TOL)

    ########################################

    # Positive scalar parameter.
    m <- Parameter(nonneg = TRUE)

    # Column vector parameter with unknown sign (by default).
    Parameter(5)

    # Matrix parameter with negative entries.
    G <- Parameter(4,7, nonpos = TRUE)

    # Assigns a constant value to G.
    value(G) <- -matrix(1, nrow = 4, ncol = 7)

    ########################################

    # Create parameter, then assign value.
    rho <- Parameter(nonneg = TRUE)
    value(rho) <- 2

    # Initialize parameter with a value.
    rho <- Parameter(nonneg = TRUE, value = 2)

    ########################################

    n <- 15
    m <- 10
    set.seed(1)
    A <- matrix(rnorm(n*m), nrow = n, ncol = m)
    b <- matrix(rnorm(n), nrow = n, ncol = 1)
    # gamma must be positive due to DCP rules.
    gamma <- Parameter(nonneg = TRUE)

    # Construct the problem.
    x <- Variable(m)
    error <- sum((A %*% x - b)^2)
    obj <- Minimize(error + gamma*p_norm(x,1))
    prob <- Problem(obj)

    # Construct a tradeoff curve of ||Ax-b||^2 vs. ||x||_1.
    sq_penalty <- c()
    l1_penalty <- c()
    x_values <- list()
    gamma_vals <- 10^seq(-4, 6, length.out = 50)
    for(val in gamma_vals) {
        value(gamma) <- val
        obj <- Minimize(error + gamma*p_norm(x,1))
        prob <- Problem(obj)
        result <- solve(prob)

        # Use value(expr) to get the numerical value of an expression in the problem.
        sq_penalty <- c(sq_penalty, result$getValue(error))
        l1_penalty <- c(l1_penalty, result$getValue(p_norm(x,1)))
        x_values <- c(x_values, list(result$getValue(x)))
    }

    ########################################
    X <- Variable(5,4)
    A <- matrix(1, nrow = 3, ncol = 5)

    # Use size(expr) to get the total number of elements.
    print(paste("size of X:", size(X)))
    print(paste("size of sum(X):", size(sum(X))))
    print(paste("size of A %*% X:", size(A %*% X)))

    # ValueError raised for invalid dimensions.
    tryCatch({
        A + X
    }, error = function(e) {
        print(e)
    })
})

test_that("Test image in-painting", {
    skip_on_cran()

    set.seed(1)
    rows <- 100
    cols <- 100

    ## Load the images and convert to arrays
    Uorig <- matrix(sample(0:255, size = rows * cols, replace = TRUE), nrow = rows, ncol = cols)

    ## Known is 1 if the pixel is known, 0 if the pixel was corrupted
    ## Known <- matrix(0, nrow = rows, ncol = cols)

    ## for(i in 1:rows) {
    ##     for(j in 1:cols) {
    ##         if(stats::runif(1) > 0.7)
    ##             Known[i,j] <- 1
    ##     }
    ## }

    Known <- sapply(seq_len(cols), function(x) as.numeric(stats::runif(n = rows) > 0.7))

    Ucorr <- Known %*% Uorig

    ## Recover the original image using total variation in-painting
    U <- Variable(rows, cols)
    obj <- Minimize(tv(U))
    constraints <- list(Known * U == Known * Ucorr)
    prob <- Problem(obj, constraints)
    res <- solve(prob, solver = "SCS")   ## Does testing on travis cause SCS to bomb?
    ## It certainly gives warnings, but that is a known problem with SCS
})

test_that("Test advanced tutorial 2", {
    skip_on_cran()

    x <- Variable()
    prob <- Problem(Minimize(x^2), list(x == 2))

    # Get ECOS arguments.
    tmp <- get_problem_data(prob, "ECOS")
    data <- tmp[[1]]
    chain <- tmp[[2]]
    inverse <- tmp[[3]]

    # Get ECOS_BB arguments.
    tmp <- get_problem_data(prob, "ECOS_BB")
    data <- tmp[[1]]
    chain <- tmp[[2]]
    inverse <- tmp[[3]]

    # Get CVXOPT arguments.
    if("CVXOPT" %in% installed_solvers()) {
        tmp <- get_problem_data(prob, "CVXOPT")
        data <- tmp[[1]]
        chain <- tmp[[2]]
        inverse <- tmp[[3]]
    }

    # Get SCS arguments.
    tmp <- get_problem_data(prob, "SCS")
    data <- tmp[[1]]
    chain <- tmp[[2]]
    inverse <- tmp[[3]]

    library(ECOSolveR)

    # Get ECOS arguments.
    tmp <- get_problem_data(prob, "ECOS")
    data <- tmp[[1]]
    chain <- tmp[[2]]
    inverse <- tmp[[3]]

    # Call ECOS solver.
    solution <- ECOSolveR::ECOS_csolve(c = data$c, G = data$G, h = data$h, dims = ECOS.dims_to_solver_dict(data$dims), A = data$A, b = data$b)

    # Unpack raw solver output.
    unpack_results(prob, solution, chain, inverse)
})

test_that("Test the log_sum_exp function", {
    skip_on_cran()

    set.seed(1)
    m <- 5
    n <- 2
    X <- matrix(1, nrow = m, ncol = n)
    w <- Variable(n)

    # expr2 <- lapply(1:m, function(i) { log_sum_exp(vstack(0, X[i,] %*% w)) })
    expr2 <- lapply(1:m, function(i) { log_sum_exp(vstack(0, matrix(X[i,], nrow = 1) %*% w)) })
    expr3 <- Reduce("+", expr2)
    obj <- Minimize(expr3)
    p <- Problem(obj)
    res <- solve(p, solver = "SCS", max_iters = 1)
})

# test_that("Test risk-return tradeoff curve", {
#     ##skip_on_cran()
#
#     library(expm)
#     n <- 4
#     S <- rbind(c( 4e-2, 6e-3,  -4e-3, 0.0),
#                c( 6e-3, 1e-2,    0.0, 0.0),
#                c(-4e-3,  0.0, 2.5e-3, 0.0),
#                c(  0.0,  0.0,    0.0, 0.0))
#     pbar <- matrix(c(0.12, 0.10, 0.07, 0.03))
#
#     N <- 100
#     Sroot <- expm::sqrtm(S)
#     x <- Variable(n, name = "x")
#     mu <- Parameter(name = "mu")
#     value(mu) <- 1   # TODO: Parameter("positive")
#     objective <- Minimize(-pbar %*% x  + mu * quad_over_lin(Sroot %*% x, 1))
#     constraints <- list(sum(x) == 1, x >= 0)
#     p <- Problem(objective, constraints)
#
#     ts <- seq(0, N-1)
#     mus <- 10^(5.0*ts/N-1.0)
#     xs <- list()
#     for(mu_val in mus) {
#         value(mu) <- mu_val
#         objective <- Minimize(-pbar %*% x  + mu * quad_over_lin(Sroot %*% x, 1))
#         p <- Problem(objective, constraints)
#         result <- solve(p)
#         xs <- c(xs, list(result$getValue(x)))
#     }
#
#     returns <- lapply(xs, function(x) { as.numeric(pbar %*% x) })
#     risks <- lapply(xs, function(x) { as.numeric(sqrt(x %*% S %*% x)) })
# })
