#ifndef STAGGEREDALGORITHM_HPP
#define STAGGEREDALGORITHM_HPP

#include "elasticity.hpp"
#include "damage.hpp"

#include "Epetra_Time.h"

namespace TPF {

  class staggeredAlgorithm : public elasticity
  {

  public:

    staggeredAlgorithm(Epetra_Comm & comm, mesh & mesh_);
    ~staggeredAlgorithm();

    void staggeredAlgorithmDirichletBC(Teuchos::ParameterList & ParametersList, bool print);

    void updateDamageHistory(Epetra_Vector & damageHistory, Epetra_Vector & displacement);

    mesh        * Mesh;
    Epetra_Comm * Comm;

    Teuchos::RCP<damage> phaseFieldBVP;

    double gc, lc;
    double lambda, mu;

    Epetra_Vector * damageSolutionOverlaped;

  };
}
#endif
