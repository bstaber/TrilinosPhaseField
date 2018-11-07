/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef BASECLASSFEM_HPP
#define BASECLASSFEM_HPP
#include "mesh.hpp"
#include "Epetra_FECrsGraph.h"
#include "Epetra_Map.h"

namespace TPF {

class baseClassFEM
{
public:
    baseClassFEM();
    ~baseClassFEM();

    mesh        * Mesh;
    Epetra_Comm * Comm;

    Epetra_Map        * OverlapMap;
    Epetra_Map        * StandardMap;
    Epetra_Import     * ImportToOverlapMap;
    Epetra_FECrsGraph * FEGraph;
};

}
#endif
