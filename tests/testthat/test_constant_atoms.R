SOLVERS_TO_TRY <- c("ECOS", "SCS")
SOLVERS_TO_TOL <- list(ECOS = 1e-6, SCS = 1e-2)

v_np <- matrix(c(-1, 2, -2), nrow = 1, ncol = 3)
LogSumExpAxis1 <- function(x) { LogSumExp(x, axis = 1) }
LogSumExpAxis2 <- function(x) { LogSumExp(x, axis = 2) }

atoms <- list(
  list(
    list(
      list(Abs, c(2,2), list(cbind(c(-5,2), c(-3,1))), Constant(cbind(c(5,2), c(3,1)))),
      list(Diag, c(2,1), list(cbind(c(-5,2), c(-3,1))), Constant(c(-5,1))),
      list(Diag, c(2,2), list(matrix(c(-5,1))), Constant(cbind(c(-5,0), c(0,1)))),
      list(Exp, c(2,2), list(cbind(c(1,0), c(2,-1))), Constant(cbind(c(exp(1),1), c(exp(2), exp(-1))))),
      list(Huber, c(2,2), list(cbind(c(0.5,-1.5), c(4,0))), Constant(cbind(c(0.25,2), c(7,0)))),
      list(function(x) { Huber(x,2.5) }, c(2,2), list(cbind(c(0.5,-1.5), c(4,0))), Constant(cbind(c(0.25,2.25), c(13.75,0)))),
      list(InvPos, c(2,2), list(cbind(c(1,2), c(3,4))), Constant(cbind(c(1,1.0/2), c(1.0/3,1.0/4)))),
      list(function(x) { (x + Constant(0))^-1 }, c(2,2), list(cbind(c(1,2), c(3,4))), Constant(cbind(c(1,1.0/2), c(1.0/3,1.0/4)))),
      list(KLDiv, c(1,1), list(exp(1), 1), Constant(1)),
      list(KLDiv, c(1,1), list(exp(1), exp(1)), Constant(0)),
      list(KLDiv, c(2,1), list(matrix(c(exp(1), 1)), 1), Constant(c(1,0))),
      list(function(x) { Kron(cbind(c(1,2), c(3,4)), x) }, c(4,4), list(cbind(c(5,6), c(7,8))),
           Constant(kronecker(cbind(c(1,2), c(3,4)), cbind(c(5,6), c(7,8))))),
      list(LambdaMax, c(1,1), list(cbind(c(2,0), c(0,1))), Constant(2)),
      list(LambdaMax, c(1,1), list(cbind(c(2,0,0), c(0,3,0), c(0,0,1))), Constant(3)),
      list(LambdaMax, c(1,1), list(cbind(c(5,7), c(7,-3))), Constant(9.06225775)),
      list(function(x) { LambdaSumLargest(x,2) }, c(1,1), list(cbind(c(1,2,3), c(2,4,5), c(3,5,6))), Constant(11.51572947)),
      list(LogSumExp, c(1,1), list(cbind(c(5,7), c(0,-3))), Constant(7.1277708268)),
      list(LogSumExpAxis1, c(2,1), list(rbind(c(5,7,1), c(0,-3,6))), Constant(c(7.12910890, 6.00259878))),
      list(LogSumExpAxis2, c(1,3), list(rbind(c(5,7,1), c(0,-3,6))), t(Constant(c(5.00671535, 7.0000454, 6.0067153)))),
      list(Logistic, c(2,2), list(cbind(c(log(5), log(7)), c(0, log(0.3)))), Constant(cbind(c(log(6),log(8)), c(log(2),log(1.3))))),
      # list(MatrixFrac, c(1,1), list(matrix(1:3), diag(3)), Constant(14)),
      # list(MatrixFrac, c(1,1), list(matrix(1:3), cbind(c(67,78,90), c(78,94,108), c(90,108,127))), Constant(0.46557377049180271)),
      # list(MatrixFrac, c(1,1), list(cbind(1:3, 4:6), cbind(c(67,78,90), c(78,94,108), c(90,108,127))), Constant(0.768852459016)),
      list(MaxElemwise, c(2,1), list(matrix(c(-5,2)), matrix(c(-3,1)), 0, matrix(c(-1,2))), Constant(c(0,2))),
      list(MaxElemwise, c(2,2), list(cbind(c(-5,2), c(-3,1)), 0, cbind(c(5,4), c(-1,2))), Constant(cbind(c(5,4), c(0,2)))),
      list(MaxEntries, c(1,1), list(cbind(c(-5,2), c(-3,1))), Constant(2)),
      list(MaxEntries, c(1,1), list(matrix(c(-5,-10))), Constant(-5)),
      list(function(x) { MaxEntries(x,axis=1) }, c(2,1), list(rbind(c(-5,2), c(-3,1))), Constant(c(2,1))),
      list(function(x) { MaxEntries(x,axis=2) }, c(1,2), list(rbind(c(-5,2), c(-3,1))), t(Constant(c(-3,2)))),
      list(function(x) { Norm(x,2) }, c(1,1), list(v_np), Constant(3)),
      list(function(x) { Norm(x,"fro") }, c(1,1), list(cbind(c(-1,2), c(3,-4))), Constant(5.47722557)),
      list(function(x) { Norm(x,1) }, c(1,1), list(v_np), Constant(5)),
      list(function(x) { Norm(x,1) }, c(1,1), list(cbind(c(-1,2), c(3,-4))), Constant(10)),
      list(function(x) { Norm(x,Inf) }, c(1,1), list(v_np), Constant(2)),
      list(function(x) { Norm(x,Inf) }, c(1,1), list(cbind(c(-1,2), c(3,-4))), Constant(4)),
      # list(function(x) { Norm(x,"nuc") }, c(1,1), list(cbind(c(2,0), c(0,1))), Constant(3)),
      # list(function(x) { Norm(x,"nuc") }, c(1,1), list(cbind(3:5, 6:8, 9:11)), Constant(23.173260452512931)),
      # list(function(x) { Norm(x,"nuc") }, c(1,1), list(cbind(3:5, 6:8)), Constant(14.618376738088918)),
      list(function(x) { SumLargest(abs(x),3) }, c(1,1), list(matrix(c(1,2,3,-4,-5))), Constant(5+4+3)),
      list(function(x) { MixedNorm(x,1,1) }, c(1,1), list(cbind(c(1,2), c(3,4), c(5,6))), Constant(21)),
      list(function(x) { MixedNorm(x,1,1) }, c(1,1), list(cbind(1:3, 4:6)), Constant(21)),
      list(function(x) { MixedNorm(x,2,1) }, c(1,1), list(cbind(c(3,3), c(4,4))), Constant(10)),
      list(function(x) { MixedNorm(x,1,Inf) }, c(1,1), list(cbind(c(1,4), c(5,6))), Constant(10)),

      list(Pnorm, c(1,1), list(matrix(1:3)), Constant(3.7416573867739413)),
      list(function(x) { Pnorm(x,1) }, c(1,1), list(matrix(c(1.1,2,-3))), Constant(6.1)),
      list(function(x) { Pnorm(x,2) }, c(1,1), list(matrix(c(1.1,2,-3))), Constant(3.7696153649941531)),
      list(function(x) { Pnorm(x,2,axis=1) }, c(2,1), list(rbind(c(1,2), c(3,4))), Constant(c(sqrt(5), 5))),
      list(function(x) { Pnorm(x,2,axis=2) }, c(1,2), list(rbind(c(1,2), c(4,5))), t(Constant(c(sqrt(17), sqrt(29))))),
      list(function(x) { Pnorm(x,Inf) }, c(1,1), list(matrix(c(1.1,2,-3))), Constant(3)),
      list(function(x) { Pnorm(x,3) }, c(1,1), list(matrix(c(1.1,2,-3))), Constant(3.3120161866074733)),
      list(function(x) { Pnorm(x,5.6) }, c(1,1), list(matrix(c(1.1,2,-3))), Constant(3.0548953718931089)),
      list(function(x) { Pnorm(x,1.2) }, c(1,1), list(cbind(1:3, 4:6)), Constant(15.971021676279573)),

      list(Pos, c(1,1), list(8), Constant(8)),
      list(Pos, c(2,1), list(matrix(c(-3,2))), Constant(c(0,2))),
      list(Neg, c(2,1), list(matrix(c(-3,3))), Constant(c(3,0))),

      list(function(x) { Power(x,0) }, c(1,1), list(7.45), Constant(1)),
      list(function(x) { Power(x,1) }, c(1,1), list(7.45), Constant(7.45)),
      list(function(x) { Power(x,2) }, c(1,1), list(7.45), Constant(55.502500000000005)),
      list(function(x) { Power(x,-1) }, c(1,1), list(7.45), Constant(0.1342281879194631)),
      list(function(x) { Power(x,-0.7) }, c(1,1), list(7.45), Constant(0.24518314363015764)),
      list(function(x) { Power(x,-1.34) }, c(1,1), list(7.45), Constant(0.06781263100321579)),
      list(function(x) { Power(x,1.34) }, c(1,1), list(7.45), Constant(14.746515290825071)),

      list(QuadOverLin, c(1,1), list(cbind(c(-1,2,-2), c(-1,2,-2)), 2), Constant(2*4.5)),
      list(QuadOverLin, c(1,1), list(v_np,2), Constant(4.5)),
      # list(function(x) { Norm(x,2) }, c(1,1), list(cbind(c(2,0), c(0,1))),  Constant(2)),
      # list(function(x) { Norm(x,2) }, c(1,1), list(cbind(3:5, 6:8, 9:11)), Constant(22.368559552680377)),
      list(function(x) { Scalene(x,2,3) }, c(2,2), list(cbind(c(-5,2), c(-3,1))), Constant(cbind(c(15,4), c(9,2)))),
      list(Square, c(2,2), list(cbind(c(-5,2), c(-3,1))), Constant(cbind(c(25,4), c(9,1)))),
      list(SumEntries, c(1,1), list(cbind(c(-5,2), c(-3,1))), Constant(-5)),
      list(function(x) { SumEntries(x,axis=2) }, c(1,2), list(cbind(c(-5,2), c(-3,1))), Constant(matrix(c(-3,-2), nrow = 1, ncol = 2))),
      list(function(x) { SumEntries(x,axis=1) }, c(2,1), list(cbind(c(-5,2), c(-3,1))), Constant(c(-8,3))),
      list(function(x) { (x + Constant(0))^2 }, c(2,2), list(cbind(c(-5,2), c(-3,1))), Constant(cbind(c(25,4), c(9,1)))),
      list(function(x) { SumLargest(x,3) }, c(1,1), list(matrix(1:5)), Constant(5+4+3)),
      list(function(x) { SumLargest(x,3) }, c(1,1), list(cbind(3:5, 6:8, 9:11)), Constant(9+10+11)),
      list(SumSquares, c(1,1), list(cbind(c(-1,2), c(3,-4))), Constant(30)),
      list(Trace, c(1,1), list(cbind(3:5, 6:8, 9:11)), Constant(3+7+11)),
      list(Trace, c(1,1), list(cbind(c(-5,2), c(-3,1))), Constant(-5+1)),
      list(TotalVariation, c(1,1), list(matrix(c(1,-1,2))), Constant(5)),
      list(TotalVariation, c(1,1), list(t(matrix(c(1,-1,2)))), Constant(5)),
      list(TotalVariation, c(1,1), list(cbind(c(-5,2), c(-3,1))), Constant(sqrt(53))),
      list(TotalVariation, c(1,1), list(cbind(c(-5,2), c(-3,1)), cbind(c(6,5), c(-4,3)), cbind(c(8,0), c(15,9))),
           Constant(base::norm(c(7,-1,-8,2,-10,7), "2"))),
      list(TotalVariation, c(1,1), list(cbind(3:5, 6:8, 9:11)), Constant(4*sqrt(10))),
      list(UpperTri, c(3,1), list(cbind(3:5, 6:8, 9:11)), Constant(c(6,9,10))),
      
      # Advanced indexing
      list(function(x) { x[cbind(c(2,3), c(1,3))] }, c(2,1), list(cbind(3:5, 6:8, 9:11)), Constant(c(4,11))),
      list(function(x) { x[c(2,3),] }, c(2,1), list(cbind(3:5, 6:8)), Constant(cbind(c(4,5), c(7,8)))),
      list(function(x) { x[cbind(3:5, 6:8) %% 2 == 0] }, c(3,1), list(cbind(3:5, 6:8)), Constant(c(4,6,8))),
      list(function(x) { x[seq(3,2,-1)] }, c(2,1), list(matrix(3:5)), Constant(c(5,4))),
      list(function(x) { x[seq(3,1,-1)] }, c(3,1), list(matrix(3:5)), Constant(c(5,4,3)))
    ),
    Minimize),
  list(
    list(
      list(Entr, c(2,2), list(cbind(c(1,exp(1)), c(exp(2), exp(-1)))), Constant(cbind(c(0,-exp(1)), c(-2*exp(2),exp(-1))))),
      list(LogDet, c(1,1), list(cbind(c(20, 8, 5, 2),
                                      c(8, 16, 2, 4),
                                      c(5, 2, 5, 2),
                                      c(2, 4, 2, 4))), Constant(7.7424020218157814)),
      # list(GeoMean, c(1,1), list(matrix(c(4,1))), Constant(2)),
      # list(GeoMean, c(1,1), list(matrix(c(0.01,7))), Constant(0.2645751311064591)),
      # list(GeoMean, c(1,1), list(matrix(c(63,7))), Constant(21)),
      # list(GeoMean, c(1,1), list(matrix(c(1,10))), Constant(sqrt(10))),
      # list(function(x) { GeoMean(x, c(1,1)) }, c(1,1), list(matrix(c(1,10))), Constant(sqrt(10))),
      # list(function(x) { GeoMean(x, c(0.4,0.8,4.9)) }, c(1,1), list(matrix(c(0.5,1.8,17))), Constant(10.04921378316062)),

      list(HarmonicMean, c(1,1), list(matrix(c(1,2,3))), Constant(1.6363636363636365)),
      list(HarmonicMean, c(1,1), list(matrix(c(2.5,2.5,2.5,2.5))), Constant(2.5)),
      list(HarmonicMean, c(1,1), list(matrix(c(0,1,2))), Constant(0)),

      list(Diff, c(2,1), list(matrix(c(1,2,3))), Constant(c(1,1))),
      list(Diff, c(1,1), list(matrix(c(1.1,2.3))), Constant(1.2)),
      list(function(x) { Diff(x, k = 2) }, c(1,1), list(matrix(c(1,2,3))), Constant(0)),
      list(Diff, c(3,1), list(matrix(c(2.1,1,4.5,-0.1))), Constant(c(-1.1,3.5,-4.6))),
      list(function(x) { Diff(x, k = 2) }, c(2,1), list(matrix(c(2.1,1,4.5,-0.1))), Constant(c(4.6,-8.1))),

      list(function(x) { Pnorm(x,0.5) }, c(1,1), list(matrix(c(1.1,2,0.1))), Constant(7.724231543909264)),
      list(function(x) { Pnorm(x,-0.4) }, c(1,1), list(matrix(c(1.1,2,0.1))), Constant(0.02713620334)),
      list(function(x) { Pnorm(x,-1) }, c(1,1), list(matrix(c(1.1,2,0.1))),  Constant(0.0876494023904)),
      list(function(x) { Pnorm(x,-2.3) }, c(1,1), list(matrix(c(1.1,2,0.1))), Constant(0.099781528576)),

      list(LambdaMin, c(1,1), list(cbind(c(2,0), c(0,1))), Constant(1)),
      list(LambdaMin, c(1,1), list(cbind(c(5,7), c(7,-3))), Constant(-7.06225775)),
      list(function(x) { LambdaSumSmallest(x,2) }, c(1,1), list(cbind(c(1,2,3), c(2,4,5), c(3,5,6))), Constant(-0.34481428)),
      list(Log, c(2,2), list(cbind(c(1,exp(1)), c(exp(2), exp(-1)))), Constant(cbind(c(0,1), c(2,-1)))),
      list(Log1p, c(2,2), list(cbind(c(0,exp(1)-1), c(exp(2)-1,exp(-1)-1))), Constant(cbind(c(0,1), c(2,-1)))),
      list(MinElemwise, c(2,1), list(matrix(c(-5,2)), matrix(c(-3,1)), 0, matrix(c(1,2))), Constant(c(-5,0))),
      list(MinElemwise, c(2,2), list(cbind(c(-5,2), c(-3,-1)), 0, cbind(c(5,4), c(-1,2))), Constant(cbind(c(-5,0), c(-3,-1)))),
      list(MinEntries, c(1,1), list(cbind(c(-5,2), c(-3,1))), Constant(-5)),
      list(MinEntries, c(1,1), list(matrix(c(-5,-10))), Constant(-10)),
      list(function(x) { x^0.25 }, c(1,1), list(7.45), Constant(7.45^0.25)),
      list(function(x) { x^0.32 }, c(2,1), list(matrix(c(7.45,3.9))), Constant(matrix(c(7.45,3.9), nrow = 2, ncol = 1)^0.32)),
      list(function(x) { x^0.9 }, c(2,2), list(cbind(c(7.45,2.2), c(4,7))), Constant(cbind(c(7.45,2.2), c(4,7))^0.9)),
      list(Sqrt, c(2,2), list(cbind(c(2,4), c(16,1))), Constant(cbind(c(1.414213562373095,2), c(4,1)))),
      list(function(x) { SumSmallest(x,3) }, c(1,1), list(matrix(c(-1,2,3,4,5))), Constant(-1+2+3)),
      list(function(x) { SumSmallest(x,4) }, c(1,1), list(cbind(c(-3,-4,5), c(6,7,8), c(9,10,11))), Constant(-3-4+5+6)),
      list(function(x) { (x + Constant(0))^0.5 }, c(2,2), list(cbind(c(2,4), c(16,1))), Constant(cbind(c(1.414213562373095,2), c(4,1))))
    ),
    Maximize
  )
)

check_solver <- function(prob, solver_name) {
  canon <- canonicalize(prob)
  objective <- canon[[1]]
  constraints <- canon[[2]]
  solver <- SOLVERS[[solver_name]]
  
  tryCatch({
      validate_solver(solver, constraints)
      return(TRUE)
    }, error = function(e) {
      return(FALSE)
  })
}

run_atom <- function(atom, problem, obj_val, solver, verbose = FALSE) {
  expect_true(is_dcp(problem))
  if(verbose) {
    print(problem@objective)
    print(problem@constraints)
    print(paste("solver", solver))
    
  }
  
  if(check_solver(problem, solver)) {
    tolerance <- SOLVERS_TO_TOL[[solver]]
    result <- solve(problem, solver = solver, verbose = verbose)

    if(tolower(result$status) %in% c("optimal", "optimal_inaccurate")) {
      if(verbose) {
        print(result$value)
        print(obj_val)
      }
      
      # TODO: Remove when Maximize inverted objective bug is fixed
      if(is(problem@objective, "Maximize"))
        diff <- (-result$value - obj_val)/(1+abs(obj_val))
      else
        diff <- (result$value - obj_val)/(1+abs(obj_val))
      expect_true(abs(diff) <= tolerance)
      
      if(abs(diff) > tolerance) {
        sink("test_constant_atoms_out.txt", append = TRUE)
        print(atom)
        cat(result$value, "\t", obj_val, "\n")
        sink()
      }
    } else
      stop("Problem status is sub-optimal: ", result$status)
  }
}

test_that("Test all constant atoms", {
  if(file.exists("test_constant_atoms_out.txt"))
    file.remove("test_constant_atoms_out.txt")
  
  for(a in atoms) {
    atom_list <- a[[1]]
    objective_type <- a[[2]]
    for(al in atom_list) {
      atom <- al[[1]]
      size <- al[[2]]
      args <- al[[3]]
      obj_val <- al[[4]]
      for(row in 1:size[1]) {
        for(col in 1:size[2]) {
          # TODO: Add more solvers as we connect them to CVXcanon
          for(solver in SOLVERS_TO_TRY) {
            # Atoms with Constant arguments
            # const_args <- lapply(args, function(arg) { Constant(arg) })
            # run_atom(atom, Problem(objective_type(do.call(atom, const_args)[row, col])),
            #         value(obj_val[row, col]), solver)
            
            # Atoms with Variable arguments
            variables <- list()
            constraints <- list()
            for(expr in args) {
              expr_size <- intf_size(expr)
              variables <- c(variables, Variable(expr_size[1], expr_size[2]))
              constraints <- c(constraints, variables[[length(variables)]] == expr)
            }
            objective <- objective_type(do.call(atom, variables)[row, col])
            print(atom)
            print(value(obj_val[row, col]))
            run_atom(atom, Problem(objective, constraints), value(obj_val[row, col]), solver)
            
            # Atoms with Parameter arguments
            # parameters <- list()
            # for(expr in args) {
            #  expr_size <- intf_size(expr)
            #  parameters <- c(parameters, Parameter(expr_size[1], expr_size[2]))
            #  value(parameters[[length(parameters)]]) <- as.matrix(expr)
            # }
            # objective <- objective_type(do.call(atom, parameters)[row, col])
            # run_atom(atom, Problem(objective), value(obj_val[row, col]), solver)
          }
        }
      }
    }
  }
})
