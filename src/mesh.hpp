#ifndef MESH_HPP
#define MESH_HPP

#include <stdio.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>
#include <iomanip>

#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_MpiComm.h"
#include "metis.h"

#include "Epetra_FECrsGraph.h"
#include "Epetra_Import.h"
#include "Epetra_Map.h"

namespace TPF {

class mesh{

  public:

    mesh(Epetra_Comm & comm, std::string & fileName_mesh, double scaling);
    ~mesh();

    int read_gmsh(std::string & filename, double scaling);
    int metis(int & NumProc);

    void partition(int & MyPID, int & NumProc);

    void get_local_nodes(int & MyPID);
    void get_cells_and_ghosts(int & MyPID);

    void shape_functions(Epetra_SerialDenseVector & N, double & xi, double & eta, double & zeta);
    void d_shape_functions(Epetra_SerialDenseMatrix & D);

    void dX_shape_functions(Epetra_SerialDenseMatrix & D, Epetra_SerialDenseMatrix JacobianMatrix, double & jac, Epetra_SerialDenseMatrix & DX);

    void jacobian_matrix(Epetra_SerialDenseMatrix & X, Epetra_SerialDenseMatrix & D, Epetra_SerialDenseMatrix & JacobianMatrix);
    void jacobian_det(Epetra_SerialDenseMatrix & JacobianMatrix, double & jac);
    void jacobian_det_faces(Epetra_SerialDenseMatrix & JacobianMatrix, double & jac);

    void set_gauss_points(int ngp, Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta);

    void computeFE();

    void build_FECrsGraphs();
    void build_EpetraMaps();

    Epetra_Comm * Comm;

    Epetra_Map        * OverlapMapU,  * OverlapMapD;
    Epetra_Map        * StandardMapU, * StandardMapD;
    Epetra_Import     * ImportToOverlapMapU, * ImportToOverlapMapD;
    Epetra_FECrsGraph * FEGraphU, * FEGraphD;

    Epetra_SerialDenseVector local_rows, vol_cells, detJac_cells;
    Epetra_SerialDenseMatrix N_cells, DX_N_cells, DY_N_cells, DZ_N_cells;

    std::vector<double> nodes_coord;
    std::vector<int>    cells_nodes;
    std::vector<int>    local_nodes_without_ghosts, local_dof_without_ghosts;
    std::vector<int>    local_nodes, local_dof;
    std::vector<int>    local_cells;

    idx_t * epart;
    idx_t * npart;
    int   * NumIndicesPerRow;

    int n_nodes = 0;
    int n_cells = 0;
    int el_type = 4;
    int n_local_nodes_without_ghosts = 0;
    int n_local_nodes = 0;
    int n_local_cells = 0;

    unsigned int n_gauss_cells;
    Epetra_SerialDenseVector gauss_weight_cells;
    Epetra_SerialDenseVector xi_cells, eta_cells, zeta_cells;
};

}
#endif
