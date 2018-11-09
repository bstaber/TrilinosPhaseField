#include "elasticity.hpp"

TPF::elasticity::elasticity(Epetra_Comm & comm, mesh & mesh_){
  Mesh = &mesh_;
  Comm = &comm;
}

TPF::elasticity::~elasticity(){

}

void TPF::elasticity::assemblePureDirichlet_homogeneousForcing(Epetra_FECrsMatrix & K){

  K.PutScalar(0.0);

  stiffness_homogeneousForcing(K);

  Comm->Barrier();

  K.GlobalAssemble();
  K.FillComplete();
}

void TPF::elasticity::stiffness_homogeneousForcing(Epetra_FECrsMatrix & K){

  int node, e_gid, error;
  int n_gauss_points = Mesh->n_gauss_cells;
  double gauss_weight;

  int *Indices_cells;
  Indices_cells = new int [3*Mesh->el_type];

  Epetra_SerialDenseMatrix Ke(3*Mesh->el_type,3*Mesh->el_type);

  Epetra_SerialDenseMatrix tangent_matrix(6,6);
  Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
  Epetra_SerialDenseMatrix matrix_B(6,3*Mesh->el_type);
  Epetra_SerialDenseMatrix B_times_TM(3*Mesh->el_type,6);

  Epetra_SerialDenseVector epsilon(6);
  Epetra_SerialDenseVector u(6);


  double phi = 0.0;

  for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
      e_gid = Mesh->local_cells[e_lid];

      for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
          node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
          dx_shape_functions(inode,0) = Mesh->DX_N_cells(inode,e_lid);
          dx_shape_functions(inode,1) = Mesh->DY_N_cells(inode,e_lid);
          dx_shape_functions(inode,2) = Mesh->DZ_N_cells(inode,e_lid);
          for (int iddl=0; iddl<3; ++iddl){
              Indices_cells[3*inode+iddl] = 3*node+iddl;
              u(3*inode+iddl) = displacementSolutionOverlaped[0][Mesh->OverlapMapU->LID(3*node+iddl)];
              for (unsigned int jnode=0; jnode<Mesh->el_type; ++jnode){
                  for (int jddl=0; jddl<3; ++jddl){
                      Ke(3*inode+iddl,3*jnode+jddl) = 0.0;
                  }
              }
          }
      }

      compute_B_matrices(dx_shape_functions, matrix_B);
      epsilon.Multiply('N', 'N', 1.0, matrix_B, u, 0.0);

      for (unsigned int gp=0; gp<n_gauss_points; ++gp){
          gauss_weight = Mesh->gauss_weight_cells(gp);

          for (unsigned inode=0; inode<4; ++inode){
            node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
            phi += Mesh->N_cells(inode,gp)*damageSolutionOverlaped[0][Mesh->OverlapMapD->LID(node)];
          }

          get_elasticity_tensor(tangent_matrix, epsilon, phi);

          error = B_times_TM.Multiply('T','N',gauss_weight*Mesh->detJac_cells(e_lid,gp),matrix_B,tangent_matrix,0.0);
          error = Ke.Multiply('N','N',1.0,B_times_TM,matrix_B,1.0);
      }

      for (unsigned int i=0; i<3*Mesh->el_type; ++i){
          for (unsigned int j=0; j<3*Mesh->el_type; ++j){
              error = K.SumIntoGlobalValues(1, &Indices_cells[i], 1, &Indices_cells[j], &Ke(i,j));
          }
      }
  }
  delete[] Indices_cells;
}

void TPF::elasticity::solve_u(Epetra_FECrsMatrix & A, Epetra_FEVector & rhs, Epetra_Vector & lhs,
                            Teuchos::ParameterList & Parameters, double & bc_disp){

  rhs.PutScalar(0.0);
  lhs.PutScalar(0.0);

  assemblePureDirichlet_homogeneousForcing(A);
  apply_dirichlet_conditions(A, rhs, bc_disp);

  int max_iter = Teuchos::getParameter<int>(Parameters.sublist("Aztec"), "AZ_max_iter");
  double tol   = Teuchos::getParameter<double>(Parameters.sublist("Aztec"), "AZ_tol");

  Epetra_LinearProblem problem(&A, &lhs, &rhs);

  AztecOO solver(problem);
  solver.SetParameters(Parameters.sublist("Aztec"));
  solver.Iterate(max_iter, tol);
}

void TPF::elasticity::compute_B_matrices(Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B){
  double factor = 1.0/std::sqrt(2.0);
  for (unsigned inode=0; inode<Mesh->el_type; ++inode){
    B(0,3*inode) = dx_shape_functions(inode,0);
    B(0,3*inode+1) = 0.0;
    B(0,3*inode+2) = 0.0;

    B(1,3*inode) = 0.0;
    B(1,3*inode+1) = dx_shape_functions(inode,1);
    B(1,3*inode+2) = 0.0;

    B(2,3*inode) = 0.0;
    B(2,3*inode+1) = 0.0;
    B(2,3*inode+2) = dx_shape_functions(inode,2);

    B(3,3*inode) = 0.0;
    B(3,3*inode+1) = factor*dx_shape_functions(inode,2);
    B(3,3*inode+2) = factor*dx_shape_functions(inode,1);

    B(4,3*inode) = factor*dx_shape_functions(inode,2);
    B(4,3*inode+1) = 0.0;
    B(4,3*inode+2) = factor*dx_shape_functions(inode,0);

    B(5,3*inode) = factor*dx_shape_functions(inode,1);
    B(5,3*inode+1) = factor*dx_shape_functions(inode,0);
    B(5,3*inode+2) = 0.0;
  }
}

int TPF::elasticity::print_solution(Epetra_Vector & solution, std::string filename){
  int NumTargetElements = 0;
  if (Comm->MyPID()==0){
      NumTargetElements = 3*Mesh->n_nodes;
  }
  Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
  Epetra_Export ExportOnRoot(*Mesh->StandardMapU,MapOnRoot);
  Epetra_MultiVector lhs_root(MapOnRoot,true);
  lhs_root.Export(solution,ExportOnRoot,Insert);
  int error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
  return error;
}
