#ifndef ECOSPROBLEM_H
#define ECOSPROBLEM_H

#include <vector>
#include <map>
#include <ecos.h>
#include "CVXcanon.hpp"
#include "LinOp.hpp"
#include "Solution.hpp"

class EcosProblem {
public:
  EcosProblem(Sense sense, LinOp *objective,
              std::map<OperatorType, std::vector<LinOp *> > constr_map,
              std::map<OperatorType, std::vector<int> > dims_map,
              std::vector<Variable> variables,
              std::map<int, int> var_offsets);

  ~EcosProblem();
  Solution solve(std::map<std::string, double> solver_options);

private:
  /* ECOS problem */
  pwork* problem;
  Sense prob_sense;
  std::vector<Variable> primal_vars;
  std::map<int, int> primal_offsets;
  std::vector<Variable> eq_dual_vars;
  std::vector<Variable> ineq_dual_vars;

  /* Dimensions */
  long n;
  long m;
  long p;
  long l;
  long ncones;
  std::vector<long> q;
  long e;

  /* Problem Matrices in CCS format */
  /* Inequality Constraints */
  std::vector<double> Gpr;
  std::vector<long> Gir;
  std::vector<long> Gjc;
  std::vector<double> h;

  /* Equality Constraints */
  std::vector<double> Apr;
  std::vector<long> Air;
  std::vector<long> Ajc;
  std::vector<double> b;

  /* Objective */
  std::vector<double> c;
  double offset;
};

#endif