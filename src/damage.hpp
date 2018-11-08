#ifndef DAMAGE_HPP
#define DAMAGE_HPP

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

#include "mesh.hpp"

namespace TPF {

  class damage
  {

  public:

    damage(Epetra_Comm & comm, mesh & mesh_, double & gc_, double & lc_);
    ~damage();

    void assemble(Epetra_FECrsMatrix & matrix, Epetra_FEVector & rhs,
                  Epetra_Vector & damageHistory);

    void solve_d(Epetra_FECrsMatrix & A, Epetra_FEVector & rhs, Epetra_Vector & lhs,
                 Teuchos::ParameterList & Parameters, Epetra_Vector & damageHistory);

    int print_solution(Epetra_Vector & lhs, std::string filename);

    mesh        * Mesh;
    Epetra_Comm * Comm;

    double gc;
    double lc;
  };

}

#endif
