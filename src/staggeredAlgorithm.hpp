#ifndef STAGGEREDALGORITHM_HPP
#define STAGGEREDALGORITHM_HPP

#include "elasticity.hpp"
#include "damage.hpp"

namespace TPF {

  class staggeredAlgorithm : public elasticity
  {

  public:

    staggeredAlgorithm(Epetra_Comm & comm, mesh & mesh_);
    ~staggeredAlgorithm();
/*
    void staggeredAlgorithmDirichletBC(Teuchos::ParameterList & ParametersList, bool print);

    void computeDisplacement(Teuchos::ParameterList & ParameterList,
                             Epetra_FECrsMatrix & matrix, Epetra_Vector & lhs, Epetra_FEVector & rhs,
                             double & bc_disp);

    void updateDamageHistory(Epetra_Vector & damageHistory,
                             Epetra_Vector & displacement);
*/

    mesh        * Mesh;
    Epetra_Comm * Comm;
    
    Teuchos::RCP<damage> phaseFieldBVP;
    Teuchos::RCP<elasticity> elasticityBVP;

    double gc, lc;
    double E, nu, lambda, mu;

    Epetra_Vector * damageSolutionOverlaped;

  };
}
#endif
