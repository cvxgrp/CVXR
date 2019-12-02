"""
Copyright 2015 Enzo Busseti
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
import scipy.sparse as sp
import mosek

def mosek_intf(A, b, G, h, c, dims, offset, solver_opts, verbose = False):
    data = {s.A : A,
            s.B : b,
            s.G : G,
            s.H : h,
            s.C : c,
            s.OFFSET : offset,
            s.DIMS : dims}

    with mosek.Env() as env:
        with env.Task(0, 0) as task:
            kwargs = sorted(solver_opts.keys())
            if "mosek_params" in kwargs:
                _handle_mosek_params(task, solver_opts["mosek_params"])
                kwargs.remove("mosek_params")
            if kwargs:
                raise ValueError("Invalid keyword-argument '%s'" % kwargs[0])
        
            if verbose:
                # Define a stream printer to grab output from MOSEK
                def streamprinter(text):
                    import sys
                    sys.stdout.write(text)
                    sys.stdout.flush()

                env.set_Stream(mosek.streamtype.log, streamprinter)
                task.set_Stream(mosek.streamtype.log, streamprinter)
            # size of problem
            numvar = len(c) + sum(dims[s.SOC_DIM])
            numcon = len(b) + dims[s.LEQ_DIM] + sum(dims[s.SOC_DIM]) + \
                     sum([el**2 for el in dims[s.SDP_DIM]])
        
            # otherwise it crashes on empty probl.
            if numvar == 0:
                result_dict = {s.STATUS: s.OPTIMAL}
                result_dict[s.PRIMAL] = []
                result_dict[s.VALUE] = 0. + data[s.OFFSET]
                result_dict[s.EQ_DUAL] = []
                result_dict[s.INEQ_DUAL] = []
                return result_dict
        
            # objective
            task.appendvars(numvar)
            task.putclist(np.arange(len(c)), c)
            task.putvarboundlist(np.arange(numvar, dtype=int),
                                 [mosek.boundkey.fr]*numvar,
                                 np.zeros(numvar),
                                 np.zeros(numvar))
            
            # SDP variables
            if sum(dims[s.SDP_DIM]) > 0:
                task.appendbarvars(dims[s.SDP_DIM])
                
            # linear equality and linear inequality constraints
            task.appendcons(numcon)
            if A.shape[0] and G.shape[0]:
                constraints_matrix = sp.bmat([[A], [G]])
            else:
                constraints_matrix = A if A.shape[0] else G
            coefficients = np.concatenate([b, h])
                    
            row, col, el = sp.find(constraints_matrix)
            task.putaijlist(row, col, el)
        
            type_constraint = [mosek.boundkey.fx] * len(b)
            type_constraint += [mosek.boundkey.up] * dims[s.LEQ_DIM]
            sdp_total_dims = sum([cdim**2 for cdim in dims[s.SDP_DIM]])
            type_constraint += [mosek.boundkey.fx] * \
                               (sum(dims[s.SOC_DIM]) + sdp_total_dims)
        
            task.putconboundlist(np.arange(numcon, dtype=int),
                                 type_constraint,
                                 coefficients,
                                 coefficients)
        
            # cones
            current_var_index = len(c)
            current_con_index = len(b) + dims[s.LEQ_DIM]
        
            for size_cone in dims[s.SOC_DIM]:
                row, col, el = sp.find(sp.eye(size_cone))
                row += current_con_index
                col += current_var_index
                task.putaijlist(row, col, el)  # add a identity for each cone
                # add a cone constraint
                task.appendcone(mosek.conetype.quad,
                                0.0,  # unused
                                np.arange(current_var_index,
                                          current_var_index + size_cone))
                current_con_index += size_cone
                current_var_index += size_cone
                
            # SDP
            for num_sdp_var, size_matrix in enumerate(dims[s.SDP_DIM]):
                for i_sdp_matrix in range(size_matrix):
                    for j_sdp_matrix in range(size_matrix):
                        coeff = 1. if i_sdp_matrix == j_sdp_matrix else .5
                        task.putbaraij(current_con_index,
                                       num_sdp_var,
                                       [task.appendsparsesymmat(size_matrix,
                                                                [max(i_sdp_matrix,
                                                                     j_sdp_matrix)],
                                                                [min(i_sdp_matrix,
                                                                     j_sdp_matrix)],
                                                                [coeff])],
                                       [1.0])
                        current_con_index += 1
        
            # solve
            task.putobjsense(mosek.objsense.minimize)
            task.optimize()
        
            if verbose:
                task.solutionsummary(mosek.streamtype.msg)

            return format_results(task, data)

def choose_solution(task):
    """Chooses between the basic and interior point solution.
    
    Parameters
    ----------
    task : mosek.Task
        The solver status interface.

    Returns
    -------
    soltype
        The preferred solution (mosek.soltype.*)
    solsta
        The status of the preferred solution (mosek.solsta.*)
    """
    import mosek
    
    def rank(status):
        # Rank solutions
        # optimal > near_optimal > anything else > None
        if status == mosek.solsta.optimal:
            return 3
        elif hasattr(mosek.solsta, 'near_optimal') and (status == mosek.solsta.near_optimal):
            return 2
        elif status is not None:
            return 1
        else:
            return 0
        
    solsta_bas, solsta_itr = None, None
        
    if task.solutiondef(mosek.soltype.bas):
        solsta_bas = task.getsolsta(mosek.soltype.bas)
        
    if task.solutiondef(mosek.soltype.itr):
        solsta_itr = task.getsolsta(mosek.soltype.itr)

    # As long as interior solution is not worse, take it
    # (for backward compatibility)
    if rank(solsta_itr) >= rank(solsta_bas):
        return mosek.soltype.itr, solsta_itr
    else:
        return mosek.soltype.bas, solsta_bas
        
def format_results(task, data):
    """Converts the solver output into standard form.

    Parameters
    ----------
    task : mosek.Task
        The solver status interface.
    data : dict
        Information about the problem.

    Returns
    -------
    dict
        The solver output in standard form.
    """

    import mosek
    # Status map is taken from:
    # https://docs.mosek.com/8.1/pythonapi/constants.html?highlight=solsta#mosek.solsta
    STATUS_MAP = {mosek.solsta.optimal: s.OPTIMAL,
                  mosek.solsta.integer_optimal: s.OPTIMAL,
                  mosek.solsta.prim_feas: s.OPTIMAL_INACCURATE,    # for integer problems
                  mosek.solsta.prim_infeas_cer: s.INFEASIBLE,
                  mosek.solsta.dual_infeas_cer: s.UNBOUNDED,
                  mosek.solsta.unknown: s.SOLVER_ERROR}
    # "Near" statuses only up to Mosek 8.1
    if hasattr(mosek.solsta, 'near_optimal'):
        STATUS_MAP_INACCURATE = {mosek.solsta.near_optimal: s.OPTIMAL_INACCURATE,
                                 mosek.solsta.near_integer_optimal: s.OPTIMAL_INACCURATE,
                                 mosek.solsta.near_prim_infeas_cer: s.INFEASIBLE_INACCURATE,
                                 mosek.solsta.near_dual_infeas_cer: s.UNBOUNDED_INACCURATE}
        STATUS_MAP.update(STATUS_MAP_INACCURATE)

    soltype, solsta = choose_solution(task)

    if solsta in STATUS_MAP:
        result_dict = {s.STATUS: STATUS_MAP[solsta]}
    else:
        result_dict = {s.STATUS: s.SOLVER_ERROR}

    # Callback data example:
    # http://docs.mosek.com/7.1/pythonapi/The_progress_call-back.html
    # Retrieving double information items:
    # http://docs.mosek.com/7.1/pythonapi/Task_getdouinf_.html#@generated-ID:5ef16e0
    # http://docs.mosek.com/7.1/pythonapi/Double_information_items.html
    result_dict[s.SOLVE_TIME] = task.getdouinf(mosek.dinfitem.optimizer_time)
    result_dict[s.SETUP_TIME] = task.getdouinf(mosek.dinfitem.presolve_time)
    result_dict[s.NUM_ITERS] = task.getintinf(mosek.iinfitem.intpnt_iter)

    if result_dict[s.STATUS] in s.SOLUTION_PRESENT:
        # get primal variables values
        result_dict[s.PRIMAL] = np.zeros(task.getnumvar(), dtype=np.float)
        task.getxx(soltype, result_dict[s.PRIMAL])
        # get obj value
        result_dict[s.VALUE] = task.getprimalobj(soltype) + \
                               data[s.OFFSET]
        # get dual
        y = np.zeros(task.getnumcon(), dtype=np.float)
        task.gety(soltype, y)
        # it appears signs are inverted
        result_dict[s.EQ_DUAL] = -y[:len(data[s.B])]
        result_dict[s.INEQ_DUAL] = \
                                   -y[len(data[s.B]):len(data[s.B])+data[s.DIMS][s.LEQ_DIM]]
        
    return result_dict

# Sets MOSEK parameters
#
# Examples of correct settings (with psolve)
#
# psolve(problem, solver = "MOSEK", mosek_params = list("MSK_IPAR_OPTIMIZER" = 1))
# psolve(problem, solver = "MOSEK", mosek_params = list("MSK_IPAR_OPTIMIZER" = "MSK_OPTIMIZER_CONIC"))
# psolve(problem, solver = "MOSEK", mosek_params = list("MSK_IPAR_OPTIMIZER" = "1"))
# psolve(problem, solver = "MOSEK", mosek_params = list("MSK_DPAR_INTPNT_CO_TOL_REL_GAP" = 1e-5))
def _handle_mosek_params(task, params):
    if params is None:
        return
    
    for param, value in params.items():
        # This is just to handle the case of integers being floats in R
        # Otherwise 1 may be formatted as 1.0 and the generic call will fail
        if '_IPAR_' in param and not '_' in str(value):
            task.putnaintparam(param, int(value))
        else:
            task.putparam(param, str(value))
