using JuMP
#using Ipopt  # Nonlinear solver
using ECOS
using Gadfly # Graphing support

N = 100   # number of chainlinks
L = 1     # difference in x-coords of endlinks
h = 2*L/N # length of each link

m = Model(solver=ECOSSolver())

@defVar(m, x[0:N])
@defVar(m, y[0:N])

# Minimize potential energy from center of mass for link
@setObjective(m, Min, sum{(y[j-1] + y[j])/2, j=1:N})

# Anchor ends
@addConstraint(m, x[0] == 0)
@addConstraint(m, y[0] == 0)
@addConstraint(m, x[N] == L)
@addConstraint(m, y[N] == 0)

# Set starting values
for j in 0:N
    setValue(x[j], 0)
    setValue(y[j], 0)
end


# Link together pieces
for j in 1:N
    @addNLConstraint(m,
        (x[j] - x[j-1])^2 + (y[j] - y[j-1])^2 <= h^2
    )
end

solve(m)

# Graph the data

x_clean = Float64[]
y_clean = Float64[]
for j in 0:N
    push!(x_clean, getValue(x)[j])
    push!(y_clean, getValue(y)[j])
end

catenary = plot(x=x_clean, y=y_clean, Coord.Cartesian(xmin=0, xmax=1,))
draw(SVG("catenary.svg", 6inch, 6inch), catenary)
