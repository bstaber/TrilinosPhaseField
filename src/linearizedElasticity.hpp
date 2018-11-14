/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef LINEARIZEDELASTICITY_HPP
#define LINEARIZEDELASTICITY_HPP

#include "linearFiniteElementProblem.hpp"

class linearizedElasticity : public linearFiniteElementProblem
{
public:
    linearizedElasticity();
    ~linearizedElasticity();

    void create_FECrsGraph();

    void aztecSolver(Epetra_FECrsMatrix & A, Epetra_FEVector & b, Epetra_Vector & u, Teuchos::ParameterList & paramList);

    void assemblePureDirichlet_homogeneousForcing(Epetra_FECrsMatrix & K);
    void stiffness_homogeneousForcing(Epetra_FECrsMatrix & K);

    void compute_B_matrices(Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B);

    int print_solution(Epetra_Vector & solution, std::string fileName);

    virtual void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, Epetra_SerialDenseMatrix & tangent_matrix) = 0;

    unsigned int n_bc_dof;
    int * dof_on_boundary;
};
#endif
