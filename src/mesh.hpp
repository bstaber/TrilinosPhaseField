#ifndef MESH_HPP
#define MESH_HPP

namespace TPF {

class mesh{

  public:

    mesh(Epetra_Comm & comm, std::string & fileName_mesh, double scaling);
    ~mesh();

    int read_gmsh(std::string & fileName_mesh, double scaling);
    int partition(int & NumProc);
    void get_local_nodes(int & MyPID);
    void get_cells_and_ghosts(int & MyPID);

    Epetra_Comm * Comm;

    std::vector<double> nodes_coord;
    std::vector<int>    cells_nodes, local_cells;
    std::vector<int>    local_nodes_without_ghosts, local_dof_without_ghosts;
    std::vector<int>    local_nodes, local_dof;

    idx_t * epart;
    idx_t * npart;
    int   * NumIndicesPerRow;

    int n_nodes = 0;
    int n_cells = 0;
    int el_type = 4;
    int n_local_nodes = 0;
    int n_local_cells = 0;
    int n_local_nodes_without_ghosts = 0;

    unsigned int n_gauss_cells;
    Epetra_SerialDenseVector gauss_weight_cells;
    Epetra_SerialDenseVector xi_cells, eta_cells, zeta_cells;

};
#endif
