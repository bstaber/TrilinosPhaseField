#ifndef ELASTICITY_HPP
#define ELASTICITY_HPP

#include "baseClassFEM.hpp"

namespace TPF {

  class elasticity
  {
  public:
    elasticity();
    ~elasticity();

    void create_FECrsGraph();

    void aztecSolver(Epetra_FECrsMatrix & A, Epetra_FEVector & b, Epetra_Vector & u, Teuchos::ParameterList & paramList);

    void assemblePureDirichlet_homogeneousForcing(Epetra_FECrsMatrix & K);
    void stiffness_homogeneousForcing(Epetra_FECrsMatrix & K);

    void compute_B_matrices(Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B);

    mesh        * Mesh;
    Epetra_Comm * Comm;

    virtual void setup_dirichlet_conditions() = 0;
    virtual void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement) = 0;
    virtual void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, Epetra_SerialDenseMatrix & tangent_matrix) = 0;

    unsigned int n_bc_dof;
    int * dof_on_boundary;
  };
}
#endif
