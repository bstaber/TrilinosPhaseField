#ifndef ELASTICITY_HPP
#define ELASTICITY_HPP

#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_VectorIn.h"

#include "Epetra_FECrsMatrix.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FEVector.h"

#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Map.h"

#include "AztecOO.h"

#include "Epetra_LAPACK.h"

#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"
#include "ml_include.h"

#include "mesh.hpp"

namespace TPF {

  class elasticity
  {
  public:
    elasticity(Epetra_Comm & comm, mesh & mesh_, double & lambda, double & mu);
    ~elasticity();

    void stiffness_homogeneousForcing(Epetra_FECrsMatrix & K, Epetra_Vector & v, Epetra_Vector & phi, Epetra_Map & OverlapMapD);

    void solve_u(Epetra_FECrsMatrix & A, Epetra_FEVector & b, Epetra_Vector & u, Epetra_Vector & v, Epetra_Vector & w, Epetra_Map & OverlapMapD, 
                 Teuchos::ParameterList & Parameters, double & bc_disp);

    void updateDamageHistory(Epetra_Vector & damageHistory, Epetra_Vector & u);

    void compute_B_matrices(Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B);

    void get_elasticity_tensor(Epetra_SerialDenseMatrix & elasticity_matrix, Epetra_SerialDenseVector & epsilon, double & phi);

    int print_solution(Epetra_Vector & solution, std::string filename);

    mesh        * Mesh;
    Epetra_Comm * Comm;

    Epetra_Map        *OverlapMap;
    Epetra_Map        *StandardMap;
    Epetra_Import     *ImportToOverlapMap;
    Epetra_FECrsGraph *FEGraph;

    double lambda = 0.0;
    double mu = 0.0;

    unsigned int n_bc_dof;
    int * dof_on_boundary;

    void setup_dirichlet_conditions();
    void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement);

    Epetra_LAPACK * Lapack;
  };
}
#endif
