/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef BASECLASSFEM_HPP
#define BASECLASSFEM_HPP

#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Epetra_FECrsMatrix.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FEVector.h"

#include "Epetra_Import.h"
#include "Epetra_Map.h"

#include "AztecOO.h"

#include "mesh.hpp"

namespace TPF {

class baseClassFEM
{
public:
    baseClassFEM(){
    };
    ~baseClassFEM(){
    };

    mesh        * Mesh;
    Epetra_Comm * Comm;

    virtual void setup_dirichlet_conditions(){
    };
    virtual void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
    };
};

}
#endif
