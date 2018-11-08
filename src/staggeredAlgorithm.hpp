#ifndef STAGGEREDALGORITHM_HPP
#define STAGGEREDALGORITHM_HPP

#include "elasticity.hpp"
#include "damage.hpp"

#include "Epetra_Time.h"
#include "Epetra_LAPACK.h"

#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"
#include "ml_include.h"

namespace TPF {

  class staggeredAlgorithm : public elasticity
  {

  public:

    staggeredAlgorithm(Epetra_Comm & comm, mesh & mesh_);
    ~staggeredAlgorithm();

    void staggeredAlgorithmDirichletBC(Teuchos::ParameterList & ParametersList, bool print);

    void updateDamageHistory(Epetra_Vector & damageHistory);

    void get_elasticity_tensor(Epetra_SerialDenseMatrix & elasticity_matrix, Epetra_SerialDenseVector & epsilon, double & phi);

    /*
    mesh        * Mesh;
    Epetra_Comm * Comm;
    */

    Teuchos::RCP<damage> phaseFieldBVP;

    double gc, lc;
    double lambda, mu;

    Epetra_Vector * damageSolutionOverlaped;
    Epetra_LAPACK * Lapack;

  };
}
#endif
