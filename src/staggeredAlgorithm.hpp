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

  class staggeredAlgorithm
  {

  public:

    staggeredAlgorithm(Epetra_Comm & comm, mesh & mesh_);
    ~staggeredAlgorithm();

    void staggeredAlgorithmDirichletBC(Teuchos::ParameterList & ParametersList, bool print);

    void updateDamageHistory(Epetra_Vector & damageHistory, Epetra_Vector & u);

    Teuchos::RCP<damage> phaseFieldBVP;

    double gc, lc;

    mesh * Mesh;
    Epetra_Comm * Comm;

  };
}
#endif
