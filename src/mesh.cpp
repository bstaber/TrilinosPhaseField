#include "mesh.hpp"

TPF::mesh::mesh(Epetra_Comm & comm, std::string & fileName_mesh, double scaling){

}

TPF::mesh::~mesh(){
  delete[] epart;
  delete[] npart;
}

int TPF::mesh::read_gmsh(std::string & filename, double scaling){

  std::ifstream meshfile;
  meshfile.open(filename.c_str());

  if (!meshfile){
      std::cout << "*ERR* Can't open the input file\n ";
      return 1;
  }

  char buf[100],c;
  double xyz[3];
  unsigned int n_total_cells;
  unsigned int el_info;
  unsigned int num;
  unsigned int nbtag;
  unsigned int tag1;
  unsigned int tag2;

  unsigned int n_cells4 = 0;
  unsigned int node;
  unsigned int nodes_cells4[4];
  std::vector<int> tetra4_nodes;

  meshfile.getline(buf,100);
  meshfile.getline(buf,100);
  meshfile.getline(buf,100);
  meshfile.getline(buf,100);
  meshfile >> n_nodes;
  meshfile.get(c);

  nodes_coord.reserve(3*n_nodes);

  for (int i=0; i<n_nodes; ++i){
      meshfile >> num;
      meshfile >> xyz[0];
      meshfile >> xyz[1];
      meshfile >> xyz[2];
      nodes_coord[3*i+0] = xyz[0]/scaling;
      nodes_coord[3*i+1] = xyz[1]/scaling;
      nodes_coord[3*i+2] = xyz[2]/scaling;
  }

  meshfile.getline(buf,100);
  meshfile.getline(buf,100);
  meshfile.getline(buf,100);
  meshfile >> n_total_cells;
  meshfile.getline(buf,100);

  for (int i=0; i<n_total_cells; ++i){
      meshfile >> num;
      meshfile >> el_info;
      meshfile >> nbtag;
      meshfile >> tag1;
      meshfile >> tag2;

      if (el_info==4){
        for (unsigned int inode=0; inode<4; ++inode){
          meshfile >> nodes_cells4[inode];
          tetra4_nodes.push_back(nodes_cells4[inode]-1);
        }
      }

  }

  meshfile.close();
  n_cells4 = tetra4_nodes.size()/4;

  if (n_cells4==0){
      std::cerr << "Your mesh is empty!\n";
  }
  else{
    n_cells = n_cells4;
    el_type = 4;
    cells_nodes.reserve(tetra4_nodes.size());
    cells_nodes = tetra4_nodes;

    gauss_points_tetra4(gauss_weight_cells, xi_cells, eta_cells, zeta_cells);
    n_gauss_cells = gauss_weight_cells.Length();
  }

  return 0;
}

int  TPF::mesh::partition(int & NumProc){
  return 0;
}

void TPF::mesh::get_local_nodes(int & MyPID){

}

void TPF::mesh::get_cells_and_ghosts(int & MyPID){

}
