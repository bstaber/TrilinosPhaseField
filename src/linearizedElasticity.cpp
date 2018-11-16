/*
Brian Staber (brian.staber@gmail.com)
*/

#include "linearizedElasticity.hpp"
#include "fepp.hpp"

linearizedElasticity::linearizedElasticity(){
}

linearizedElasticity::~linearizedElasticity(){
  delete[] dof_on_boundary;
}

void linearizedElasticity::create_FECrsGraph(){

  FEGraph = new Epetra_FECrsGraph(Copy,*StandardMap,100);
  int eglob, node;
  int *index;
  index = new int [3*Mesh->el_type];

  for (int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
      eglob = Mesh->local_cells[e_lid];
      for (int inode=0; inode<Mesh->el_type; ++inode){
          node = Mesh->cells_nodes[Mesh->el_type*eglob+inode];
          for (int ddl=0; ddl<3; ++ddl){
              index[3*inode+ddl] = 3*node+ddl;
          }
      }
      for (int i=0; i<3*Mesh->el_type; ++i){
          for (int j=0; j<3*Mesh->el_type; ++j){
              FEGraph->InsertGlobalIndices(1, &index[i], 1, &index[j]);
          }
      }

  }
  Comm->Barrier();
  FEGraph->GlobalAssemble();
  delete [] index;
}

void linearizedElasticity::aztecSolver(Epetra_FECrsMatrix & A, Epetra_FEVector & b, Epetra_Vector & u, Teuchos::ParameterList & paramList){
  u.PutScalar(0.0);
  Epetra_LinearProblem problem;
  AztecOO solver;
  problem.SetOperator(&A);
  problem.SetLHS(&u);
  problem.SetRHS(&b);
  solver.SetProblem(problem);
  solver.SetParameters(paramList);
  double tol  = Teuchos::getParameter<double>(paramList,"AZ_tol");
  int maxIter = Teuchos::getParameter<int>(paramList,"AZ_max_iter");
  solver.Iterate(maxIter,tol);
}

void linearizedElasticity::assemblePureDirichlet_homogeneousForcing(Epetra_FECrsMatrix & K, Epetra_Vector & phi, Epetra_Map & OverlapMapD){

  int error;

  K.PutScalar(0.0);

  stiffness_homogeneousForcing(K, phi, OverlapMapD);

  Comm->Barrier();

  error=K.GlobalAssemble();
  error=K.FillComplete();
}

void linearizedElasticity::stiffness_homogeneousForcing(Epetra_FECrsMatrix & K, Epetra_Vector & phi, Epetra_Map & OverlapMapD){

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
  Epetra_SerialDenseVector epsilon(6), cell_u(6);
  Epetra_SerialDenseVector shape_functions(Mesh->el_type);

  double xi, eta, zeta;

  double phi_e;
  for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
      e_gid = Mesh->local_cells[e_lid];

      for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
          node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
          for (int iddl=0; iddl<3; ++iddl){
              Indices_cells[3*inode+iddl] = 3*node+iddl;
              for (unsigned int jnode=0; jnode<Mesh->el_type; ++jnode){
                  for (int jddl=0; jddl<3; ++jddl){
                      Ke(3*inode+iddl,3*jnode+jddl) = 0.0;
                  }
              }
          }
      }

      for (unsigned int gp=0; gp<n_gauss_points; ++gp){
          gauss_weight = Mesh->gauss_weight_cells(gp);
          xi   = Mesh->xi_cells[gp];
          eta  = Mesh->eta_cells[gp];
          zeta = Mesh->zeta_cells[gp];
          tetra4::shape_functions(shape_functions, xi, eta, zeta);
          phi_e = 0.0;
          for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
              dx_shape_functions(inode,0) = Mesh->DX_N_cells(gp+n_gauss_points*inode,e_lid);
              dx_shape_functions(inode,1) = Mesh->DY_N_cells(gp+n_gauss_points*inode,e_lid);
              dx_shape_functions(inode,2) = Mesh->DZ_N_cells(gp+n_gauss_points*inode,e_lid);

              node   = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
              phi_e += shape_functions(inode)*phi[OverlapMapD.LID(node)];
          }
          compute_B_matrices(dx_shape_functions,matrix_B);
          get_elasticity_tensor(e_lid, gp, phi_e, epsilon, tangent_matrix);

          error = B_times_TM.Multiply('T','N',gauss_weight*Mesh->detJac_cells(e_lid,gp),matrix_B,tangent_matrix,0.0);
          error = Ke.Multiply('N','N',1.0,B_times_TM,matrix_B,1.0);
      }

      for (unsigned int i=0; i<3*Mesh->el_type; ++i){
          for (unsigned int j=0; j<3*Mesh->el_type; ++j){
              error = K.SumIntoGlobalValues(1, &Indices_cells[i], 1, &Indices_cells[j], &Ke(i,j));
          }
      }
  }
  delete [] Indices_cells;
}

void linearizedElasticity::compute_B_matrices(Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B){
  double factor = 1.0; ///std::sqrt(2.0);
  for (unsigned inode=0; inode<Mesh->el_type; ++inode){
      B(0,3*inode)   = dx_shape_functions(inode,0);
      B(0,3*inode+1) = 0.0;
      B(0,3*inode+2) = 0.0;

      B(1,3*inode)   = 0.0;
      B(1,3*inode+1) = dx_shape_functions(inode,1);
      B(1,3*inode+2) = 0.0;

      B(2,3*inode)   = 0.0;
      B(2,3*inode+1) = 0.0;
      B(2,3*inode+2) = dx_shape_functions(inode,2);

      B(3,3*inode)   = 0.0;
      B(3,3*inode+1) = factor*dx_shape_functions(inode,2);
      B(3,3*inode+2) = factor*dx_shape_functions(inode,1);

      B(4,3*inode)   = factor*dx_shape_functions(inode,2);
      B(4,3*inode+1) = 0.0;
      B(4,3*inode+2) = factor*dx_shape_functions(inode,0);

      B(5,3*inode)   = factor*dx_shape_functions(inode,1);
      B(5,3*inode+1) = factor*dx_shape_functions(inode,0);
      B(5,3*inode+2) = 0.0;
  }
}

int linearizedElasticity::print_solution(Epetra_Vector & solution, std::string fileName){

  int NumTargetElements = 0;
  if (Comm->MyPID()==0){
      NumTargetElements = 3*Mesh->n_nodes;
  }
  Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
  Epetra_Export ExportOnRoot(*StandardMap,MapOnRoot);
  Epetra_MultiVector lhs_root(MapOnRoot,true);
  lhs_root.Export(solution,ExportOnRoot,Insert);

  int error = EpetraExt::MultiVectorToMatrixMarketFile(fileName.c_str(),lhs_root,0,0,false);

  return error;
}
