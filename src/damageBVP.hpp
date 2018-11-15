/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef DAMAGEBVP_HPP
#define DAMAGEBVP_HPP

#include "linearFiniteElementProblem.hpp"

class damageBVP : public linearFiniteElementProblem
{
public:
    double gc;
    double lc;

    damageBVP(Epetra_Comm & comm, mesh & mesh, double & gc_, double & lc_);
    ~damageBVP();

    void solve(Teuchos::ParameterList & Parameters,
               Epetra_FECrsMatrix & matrix, Epetra_Vector & lhs, Epetra_FEVector & rhs,
               Epetra_Vector & damageHistory, Epetra_Map & GaussMap);

    void assemble(Epetra_FECrsMatrix & matrix, Epetra_FEVector & rhs,
                  Epetra_Vector & damageHistory, Epetra_Map & GaussMap);

    void create_FECrsGraph();
    int print_solution(Epetra_Vector & lhs, std::string fileName);
    void setup_dirichlet_conditions();
    void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement);
};

#endif
