> x <- Variable(1)
> obj <- Maximize(Log(x))
> prob <- Problem(obj)
> get_problem_data(prob, "ECOS")
Before Build Matrix 1
Processing constraint 0
After constructing external ptr
Instantiating ProblemData-R6
Before Build Matrix 1
Processing constraint 0
After constructing external ptr
Instantiating ProblemData-R6
$c
[1]  0 -1

$offset
[1] 0

$A
0 x 2 sparse Matrix of class "ngCMatrix"


$b
numeric(0)

$G
3 x 2 sparse Matrix of class "dgCMatrix"

[1,]  . -1
[2,]  .  .
[3,] -1  .

$h
[1] 0 1 0

$dims
$dims$f
[1] 0

$dims$l
[1] 0

$dims$ep
[1] 1


$bool_vars_idx
list()

$int_vars_idx
list()

>
