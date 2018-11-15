/*
Brian Staber (brian.staber@gmail.com)
*/

#include "elasticBVP.hpp"

elasticBVP::elasticBVP(Epetra_Comm & comm, mesh & mesh_, Teuchos::ParameterList & Parameters){

  Comm = &comm;
  Mesh = &mesh_;

  StandardMap = new Epetra_Map(-1, 3*Mesh->n_local_nodes_without_ghosts, &Mesh->local_dof_without_ghosts[0], 0, *Comm);
  OverlapMap  = new Epetra_Map(-1, 3*Mesh->n_local_nodes, &Mesh->local_dof[0], 0, *Comm);
  ImportToOverlapMap = new Epetra_Import(*OverlapMap, *StandardMap);
  create_FECrsGraph();

}

elasticBVP::~elasticBVP(){
}

void elasticBVP::computeDisplacement(Teuchos::ParameterList & Parameters,
                                     Epetra_FECrsMatrix & matrix, Epetra_Vector & lhs, Epetra_FEVector & rhs,
                                     Epetra_Vector & phi, Epetra_Map & OverlapMapD, double & bc_disp){

  rhs.PutScalar(0.0);
  lhs.PutScalar(0.0);

  assemblePureDirichlet_homogeneousForcing(matrix, phi, OverlapMapD);
  apply_dirichlet_conditions(matrix, rhs, bc_disp);

  int max_iter = Teuchos::getParameter<int>(Parameters.sublist("Aztec"), "AZ_max_iter");
  double tol   = Teuchos::getParameter<double>(Parameters.sublist("Aztec"), "AZ_tol");

  Epetra_LinearProblem problem(&matrix, &lhs, &rhs);

  AztecOO solver(problem);
  solver.SetParameters(Parameters.sublist("Aztec"));
  solver.Iterate(max_iter, tol);
}

Epetra_Map elasticBVP::constructGaussMap(){
  int e_gid;
  int n_local_cells = Mesh->n_local_cells;
  int n_gauss_cells = Mesh->n_gauss_cells;
  std::vector<int> local_gauss_points(n_local_cells*n_gauss_cells);
  for (unsigned int e_lid=0; e_lid<n_local_cells; ++e_lid){
      e_gid = Mesh->local_cells[e_lid];
      for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
          local_gauss_points[e_lid*n_gauss_cells+gp] = e_gid*n_gauss_cells+gp;
      }

  }
  Epetra_Map GaussMap(-1, n_local_cells*n_gauss_cells, &local_gauss_points[0], 0, *Comm);
  return GaussMap;
}
