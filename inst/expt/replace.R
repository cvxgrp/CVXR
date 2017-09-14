library(cvxr)
x <- Variable()
obj <- Minimize(Log(x))
constr <- list(x >= 1)
prob <- Problem(obj, constr)

value(x)    # This will be NA

value(variables(prob@objective)) <- 5



                                        # Set value of all variables to 5
for(var_ in variables(prob))
    value(var_) <- 5

# This should print 5, but still gives me NA
for(var_ in variables(prob))
    print(value(var_))

vars <- variables(prob@objective)
constrs_ <- lapply(object@constraints, function(constr) { variables(constr) })
unique(flatten_list(c(vars_, constrs_)))   # Remove duplicates
