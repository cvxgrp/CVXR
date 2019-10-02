## These tests were moved here because scs causes errors on R-devel on travis
## Also because they go over the time limit
TOL <- 1e-6

test_that("Test Chebyshev center", {
    # Find the largest Euclidean ball in the polyhedron.
    # Goal is to find the largest Euclidean ball (i.e. its center and radius)
    # that lies in a polyhedron described by linear inequalities in this fashion:
    # P = {x : a_i'*x <= b_i, i=1,...,m} where x is in R^2.
    skip_on_cran()
    
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
    ## Problem data
    skip_on_cran()
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
    res <- solve(prob, solver = "SCS") ## Does testing on travis cause SCS to bomb?
    ## It certainly gives warnings, but that is a known problem with SCS
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
