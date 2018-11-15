/*
Brian Staber (brian.staber@gmail.com)
*/

#include "elasticBVP.hpp"

elasticBVP::elasticBVP(Epetra_Comm & comm, mesh & mesh_, Teuchos::ParameterList & Parameters){

  Comm = &comm;
  Mesh = &mesh_;

  gc = Teuchos::getParameter<double>(Parameters.sublist("Damage"), "gc");
  lc = Teuchos::getParameter<double>(Parameters.sublist("Damage"), "lc");
  E  = Teuchos::getParameter<double>(Parameters.sublist("Elasticity"), "young");
  nu = Teuchos::getParameter<double>(Parameters.sublist("Elasticity"), "poisson");

  damageInterface = Teuchos::rcp(new damageBVP(comm, *Mesh, gc, lc));

  StandardMap = new Epetra_Map(-1, 3*Mesh->n_local_nodes_without_ghosts, &Mesh->local_dof_without_ghosts[0], 0, *Comm);
  OverlapMap  = new Epetra_Map(-1, 3*Mesh->n_local_nodes, &Mesh->local_dof[0], 0, *Comm);
  ImportToOverlapMap = new Epetra_Import(*OverlapMap, *StandardMap);
  create_FECrsGraph();

  lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
  mu     = E/(2.0*(1.0+nu));

  elasticity.Reshape(6,6);
  double c11 = E*(1.0-nu)/((1.0+nu)*(1.0-2.0*nu));
  double c12 = E*nu/((1.0+nu)*(1.0-2.0*nu));
  double c44 = E/(2.0*(1.0+nu));

  elasticity(0,0) = c11; elasticity(0,1) = c12; elasticity(0,2) = c12; elasticity(0,3) = 0.0; elasticity(0,4) = 0.0; elasticity(0,5) = 0.0;
  elasticity(1,0) = c12; elasticity(1,1) = c11; elasticity(1,2) = c12; elasticity(1,3) = 0.0; elasticity(1,4) = 0.0; elasticity(1,5) = 0.0;
  elasticity(2,0) = c12; elasticity(2,1) = c12; elasticity(2,2) = c11; elasticity(2,3) = 0.0; elasticity(2,4) = 0.0; elasticity(2,5) = 0.0;
  elasticity(3,0) = 0.0; elasticity(3,1) = 0.0; elasticity(3,2) = 0.0; elasticity(3,3) = c44; elasticity(3,4) = 0.0; elasticity(3,5) = 0.0;
  elasticity(4,0) = 0.0; elasticity(4,1) = 0.0; elasticity(4,2) = 0.0; elasticity(4,3) = 0.0; elasticity(4,4) = c44; elasticity(4,5) = 0.0;
  elasticity(5,0) = 0.0; elasticity(5,1) = 0.0; elasticity(5,2) = 0.0; elasticity(5,3) = 0.0; elasticity(5,4) = 0.0; elasticity(5,5) = c44;

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

void elasticBVP::updateDamageHistory(Epetra_Vector & damageHistory,
                                     Epetra_Vector & displacement,
                                     Epetra_Map & GaussMap){

  Epetra_Vector u(*OverlapMap);
  u.Import(displacement, *ImportToOverlapMap, Insert);

  int n_gauss_points = Mesh->n_gauss_cells;

  double trepsilon, trepsilon2, potential, history;

  Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
  Epetra_SerialDenseMatrix matrix_B(6,3*Mesh->el_type);
  Epetra_SerialDenseVector cells_u(3*Mesh->el_type);
  Epetra_SerialDenseVector epsilon(6);

  int egid, node, id;
  for (unsigned int elid=0; elid<Mesh->n_local_cells; ++elid){
    egid = Mesh->local_cells[elid];
    for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
      node = Mesh->cells_nodes[Mesh->el_type*egid+inode];
      for (unsigned int ddl=0; ddl<3; ++ddl){
        id = 3*node+ddl;
        cells_u(3*inode+ddl) = u[OverlapMap->LID(id)];
      }
    }
    for (unsigned int gp=0; gp<n_gauss_points; ++gp){
        for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
            dx_shape_functions(inode,0) = Mesh->DX_N_cells(gp+n_gauss_points*inode,elid);
            dx_shape_functions(inode,1) = Mesh->DY_N_cells(gp+n_gauss_points*inode,elid);
            dx_shape_functions(inode,2) = Mesh->DZ_N_cells(gp+n_gauss_points*inode,elid);
        }

        compute_B_matrices(dx_shape_functions,matrix_B);
        epsilon.Multiply('N','N',1.0,matrix_B,cells_u,0.0);

        trepsilon  = epsilon(0) + epsilon(1) + epsilon(2);
        trepsilon2 = epsilon(0)*epsilon(0) + epsilon(1)*epsilon(1) + epsilon(2)*epsilon(2) +
                     0.5*epsilon(3)*epsilon(3) + 0.5*epsilon(4)*epsilon(4) + 0.5*epsilon(5)*epsilon(5);

        id = n_gauss_points*egid+gp;
        history   = damageHistory[GaussMap.LID(id)];
        potential = (lambda/2.0)*trepsilon*trepsilon + mu*trepsilon2;
        if (potential>history){
          damageHistory[GaussMap.LID(id)] = potential;
        }
    }
  }

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
