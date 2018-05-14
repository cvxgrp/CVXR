"""
Copyright 2017 Steven Diamond
 [ Modified by Balasubramanian Narasimhan ]

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import settings as s
import numpy as np
from scipy.sparse import dok_matrix

def gurobi_intf(Amat, b, c, bool_idx, int_idx, dims, offset, solver_opts, verbose = False):
    import gurobipy

    STATUS_MAP = {2: s.OPTIMAL,
                  3: s.INFEASIBLE,
                  5: s.UNBOUNDED,
                  4: s.SOLVER_ERROR,
                  6: s.SOLVER_ERROR,
                  7: s.SOLVER_ERROR,
                  8: s.SOLVER_ERROR,
                  # TODO could be anything.
                  # means time expired.
                  9: s.OPTIMAL_INACCURATE,
                  10: s.SOLVER_ERROR,
                  11: s.SOLVER_ERROR,
                  12: s.SOLVER_ERROR,
                  13: s.SOLVER_ERROR}
    
    A = dok_matrix(Amat)
    # Save the dok_matrix.
    data = {s.A : A,
            s.B : b,
            s.C : c,
            s.BOOL_IDX : bool_idx,
            s.INT_IDX : int_idx,
            s.OFFSET : offset,
            s.DIMS : dims}

    n = c.shape[0]

    model = gurobipy.Model()
    variables = []
    for i in range(n):
        # Set variable type.
        if i in data[s.BOOL_IDX]:
            vtype = gurobipy.GRB.BINARY
        elif i in data[s.INT_IDX]:
            vtype = gurobipy.GRB.INTEGER
        else:
            vtype = gurobipy.GRB.CONTINUOUS
        variables.append(
            model.addVar(
                obj=c[i],
                name="x_%d" % i,
                vtype=vtype,
                # Gurobi's default LB is 0 (WHY???)
                lb=-gurobipy.GRB.INFINITY,
                ub=gurobipy.GRB.INFINITY)
        )
    model.update()

    eq_constrs = add_model_lin_constr(model, variables,
                                      range(data[s.DIMS][s.EQ_DIM]),
                                      gurobipy.GRB.EQUAL,
                                      A, b)
    leq_start = data[s.DIMS][s.EQ_DIM]
    leq_end = data[s.DIMS][s.EQ_DIM] + data[s.DIMS][s.LEQ_DIM]
    ineq_constrs = add_model_lin_constr(model, variables,
                                        range(leq_start, leq_end),
                                        gurobipy.GRB.LESS_EQUAL,
                                        A, b)
    soc_start = leq_end
    soc_constrs = []
    new_leq_constrs = []
    for constr_len in data[s.DIMS][s.SOC_DIM]:
        soc_end = soc_start + constr_len
        soc_constr, new_leq, new_vars = add_model_soc_constr(
            model, variables, range(soc_start, soc_end),
            A, b
        )
        soc_constrs.append(soc_constr)
        new_leq_constrs += new_leq
        variables += new_vars
        soc_start += constr_len
        
    gur_constrs = eq_constrs + ineq_constrs + \
                  soc_constrs + new_leq_constrs
    model.update()

    # Set verbosity and other parameters
    model.setParam("OutputFlag", verbose)
    # TODO user option to not compute duals.
    model.setParam("QCPDual", True)

    for key, value in solver_opts.items():
        model.setParam(key, value)

    results_dict = {}
    try:
        model.optimize()
        results_dict["primal objective"] = model.ObjVal
        results_dict["x"] = np.array([v.X for v in variables])

        # Only add duals if not a MIP.
        # Not sure why we need to negate the following,
        # but need to in order to be consistent with other solvers.
        if not is_mip(data):
            vals = []
            for lc in gur_constrs:
                if lc is not None:
                    if isinstance(lc, gurobipy.QConstr):
                        vals.append(lc.QCPi)
                    else:
                        vals.append(lc.Pi)
                else:
                    vals.append(0)
            results_dict["y"] = -np.array(vals)

        results_dict["status"] = STATUS_MAP.get(model.Status, s.SOLVER_ERROR)
    except gurobipy.GurobiError:
        results_dict["status"] = s.SOLVER_ERROR

    results_dict["model"] = model
    results_dict["variables"] = variables
    results_dict["gur_constrs"] = gur_constrs
    results_dict[s.SOLVE_TIME] = model.Runtime

    return format_results(results_dict, data)

def add_model_lin_constr(model, variables,
                         rows, ctype,
                         mat, vec):
    """Adds EQ/LEQ constraints to the model using the data from mat and vec.

    Parameters
    ----------
    model : GUROBI model
        The problem model.
    variables : list
        The problem variables.
    rows : range
        The rows to be constrained.
    ctype : GUROBI constraint type
        The type of constraint.
    mat : SciPy COO matrix
        The matrix representing the constraints.
    vec : NDArray
        The constant part of the constraints.

    Returns
    -------
    list
        A list of constraints.
    """
    import gurobipy
    # print(variables)
    # print(rows)
    # print(ctype)
    # print(mat)
    # print(vec)
    constr = []
    expr_list = {i: [] for i in rows}
    for (i, j), c in mat.iteritems():
        v = variables[j]
        if i in rows:
            expr_list[i].append((c, v))

    for i in rows:
        # Ignore empty constraints.
        if expr_list[i]:
            expr = gurobipy.LinExpr(expr_list[i])
            constr.append(
                model.addConstr(expr, ctype, vec[i])
            )
        else:
            constr.append(None)
    return constr

def add_model_soc_constr(model, variables,
                         rows, mat, vec):
    """Adds SOC constraint to the model using the data from mat and vec.

    Parameters
    ----------
    model : GUROBI model
        The problem model.
    variables : list
        The problem variables.
    rows : range
        The rows to be constrained.
    mat : SciPy COO matrix
        The matrix representing the constraints.
    vec : NDArray
        The constant part of the constraints.

    Returns
    -------
    tuple
        A tuple of (QConstr, list of Constr, and list of variables).
    """
    import gurobipy
    # Assume first expression (i.e. t) is nonzero.
    expr_list = {i: [] for i in rows}
    for (i, j), c in mat.iteritems():
        v = variables[j]
        if i in rows:
            expr_list[i].append((c, v))

    lin_expr_list = [vec[i] - gurobipy.LinExpr(expr_list[i]) for i in rows]

    # Make a variable and equality constraint for each term.
    soc_vars = [
        model.addVar(
            obj=0,
            name="soc_t_%d" % rows[0],
            vtype=gurobipy.GRB.CONTINUOUS,
            lb=0,
            ub=gurobipy.GRB.INFINITY)
    ]
    for i in rows[1:]:
        soc_vars += [
            model.addVar(
                obj=0,
                name="soc_x_%d" % i,
                vtype=gurobipy.GRB.CONTINUOUS,
                lb=-gurobipy.GRB.INFINITY,
                ub=gurobipy.GRB.INFINITY)
        ]
    model.update()

    new_lin_constrs = []
    for i, _ in enumerate(lin_expr_list):
        new_lin_constrs += [
            model.addConstr(soc_vars[i] == lin_expr_list[i])
        ]

    t_term = soc_vars[0]*soc_vars[0]
    x_term = gurobipy.quicksum([var*var for var in soc_vars[1:]])
    return (model.addQConstr(x_term <= t_term),
            new_lin_constrs,
            soc_vars)

def format_results(results_dict, data):
    """Converts the solver output into standard form.

    Parameters
    ----------
    results_dict : dict
        The solver output.
    data : dict
        Information about the problem.
    cached_data : dict
        A map of solver name to cached problem data.

    Returns
    -------
    dict
        The solver output in standard form.
    """
    dims = data[s.DIMS]
    # if results_dict["status"] != s.SOLVER_ERROR:
    #     solver_cache = cached_data[self.name()]
    #     solver_cache.prev_result = {
    #         "model": results_dict["model"],
    #         "variables": results_dict["variables"],
    #         "gur_constrs": results_dict["gur_constrs"],
    #         "c": data[s.C],
    #         "A": data[s.A],
    #         "b": data[s.B],
    #     }
    new_results = {}
    new_results[s.STATUS] = results_dict['status']
    new_results[s.SOLVE_TIME] = results_dict[s.SOLVE_TIME]
    if new_results[s.STATUS] in s.SOLUTION_PRESENT:
        primal_val = results_dict['primal objective']
        new_results[s.VALUE] = primal_val + data[s.OFFSET]
        new_results[s.PRIMAL] = results_dict['x']
        if not is_mip(data):
            new_results[s.EQ_DUAL] = results_dict["y"][0:dims[s.EQ_DIM]]
            new_results[s.INEQ_DUAL] = results_dict["y"][dims[s.EQ_DIM]:]

    return new_results

def is_mip(data):
    """Is the problem a mixed integer program?
    """
    return len(data[s.BOOL_IDX]) > 0 or len(data[s.INT_IDX]) > 0
