#include "elasticity.hpp"

TPF::elasticity::elasticity(Epetra_Comm & comm, mesh & mesh_, double & lambda_, double & mu_){
  Mesh = &mesh_;
  Comm = &comm;

  Lapack = new Epetra_LAPACK();

  lambda = lambda_;
  mu = mu_;

  setup_dirichlet_conditions();
}

TPF::elasticity::~elasticity(){

}

void TPF::elasticity::stiffness_homogeneousForcing(Epetra_FECrsMatrix & K, Epetra_Vector & v, Epetra_Vector & phi){

  K.PutScalar(0.0);

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
  Epetra_SerialDenseVector u(3*Mesh->el_type);

  double pf;
  for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
      e_gid = Mesh->local_cells[e_lid];

      for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
          node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
          dx_shape_functions(inode,0) = Mesh->DX_N_cells(inode,e_lid);
          dx_shape_functions(inode,1) = Mesh->DY_N_cells(inode,e_lid);
          dx_shape_functions(inode,2) = Mesh->DZ_N_cells(inode,e_lid);

          for (int iddl=0; iddl<3; ++iddl){
              Indices_cells[3*inode+iddl] = 3*node+iddl;
              u(3*inode+iddl) = v[Mesh->OverlapMapU->LID(3*node+iddl)];
              for (unsigned int jnode=0; jnode<Mesh->el_type; ++jnode){
                  for (int jddl=0; jddl<3; ++jddl){
                      Ke(3*inode+iddl,3*jnode+jddl) = 0.0;
                  }
              }
          }

      }

      compute_B_matrices(dx_shape_functions, matrix_B);
      epsilon.Multiply('N', 'N', 1.0, matrix_B, u, 0.0);

      pf = 0.0;
      for (unsigned int gp=0; gp<n_gauss_points; ++gp){
          gauss_weight = Mesh->gauss_weight_cells(gp);

          for (unsigned inode=0; inode<4; ++inode){
            node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
            pf  += Mesh->N_cells(inode,gp)*phi[Mesh->OverlapMapD->LID(node)];
          }

          get_elasticity_tensor(tangent_matrix, epsilon, pf);

          error = B_times_TM.Multiply('T', 'N', gauss_weight*Mesh->detJac_cells(e_lid), matrix_B, tangent_matrix, 0.0);
          error = Ke.Multiply('N', 'N', 1.0, B_times_TM, matrix_B, 1.0);
      }

      for (unsigned int i=0; i<3*Mesh->el_type; ++i){
          for (unsigned int j=0; j<3*Mesh->el_type; ++j){
              error = K.SumIntoGlobalValues(1, &Indices_cells[i], 1, &Indices_cells[j], &Ke(i,j));
          }
      }

  }

  delete[] Indices_cells;

  Comm->Barrier();

  K.GlobalAssemble();
  K.FillComplete();
}

void TPF::elasticity::solve_u(Epetra_FECrsMatrix & A, Epetra_FEVector & rhs, Epetra_Vector & lhs, Epetra_Vector & v, Epetra_Vector & w,
                              Teuchos::ParameterList & Parameters, double & bc_disp){

  stiffness_homogeneousForcing(A, v, w);
  apply_dirichlet_conditions(A, rhs, bc_disp);

  int max_iter = Teuchos::getParameter<int>(Parameters.sublist("Elasticity").sublist("Aztec"), "AZ_max_iter");
  double tol   = Teuchos::getParameter<double>(Parameters.sublist("Elasticity").sublist("Aztec"), "AZ_tol");

  Epetra_LinearProblem problem(&A, &lhs, &rhs);
  AztecOO solver(problem);

  solver.SetParameters(Parameters.sublist("Elasticity").sublist("Aztec"));
  solver.Iterate(max_iter, tol);
}

void TPF::elasticity::compute_B_matrices(Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B){

  for (unsigned inode=0; inode<Mesh->el_type; ++inode){
    B(0,3*inode+0) = dx_shape_functions(inode,0);
    B(0,3*inode+1) = 0.0;
    B(0,3*inode+2) = 0.0;

    B(1,3*inode+0) = 0.0;
    B(1,3*inode+1) = dx_shape_functions(inode,1);
    B(1,3*inode+2) = 0.0;

    B(2,3*inode+0) = 0.0;
    B(2,3*inode+1) = 0.0;
    B(2,3*inode+2) = dx_shape_functions(inode,2);

    B(3,3*inode+0) = 0.0;
    B(3,3*inode+1) = dx_shape_functions(inode,2);
    B(3,3*inode+2) = dx_shape_functions(inode,1);

    B(4,3*inode+0) = dx_shape_functions(inode,2);
    B(4,3*inode+1) = 0.0;
    B(4,3*inode+2) = dx_shape_functions(inode,0);

    B(5,3*inode+0) = dx_shape_functions(inode,1);
    B(5,3*inode+1) = dx_shape_functions(inode,0);
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

void TPF::elasticity::get_elasticity_tensor(Epetra_SerialDenseMatrix & elasticity_matrix, Epetra_SerialDenseVector & epsilon, double & phi){

  double c11 = lambda+2.0*mu;
  double c12 = lambda;
  double c44 = mu;

  elasticity_matrix(0,0) = c11; elasticity_matrix(0,1) = c12; elasticity_matrix(0,2) = c12; elasticity_matrix(0,3) = 0.0; elasticity_matrix(0,4) = 0.0; elasticity_matrix(0,5) = 0.0;
  elasticity_matrix(1,0) = c12; elasticity_matrix(1,1) = c11; elasticity_matrix(1,2) = c12; elasticity_matrix(1,3) = 0.0; elasticity_matrix(1,4) = 0.0; elasticity_matrix(1,5) = 0.0;
  elasticity_matrix(2,0) = c12; elasticity_matrix(2,1) = c12; elasticity_matrix(2,2) = c11; elasticity_matrix(2,3) = 0.0; elasticity_matrix(2,4) = 0.0; elasticity_matrix(2,5) = 0.0;
  elasticity_matrix(3,0) = 0.0; elasticity_matrix(3,1) = 0.0; elasticity_matrix(3,2) = 0.0; elasticity_matrix(3,3) = c44; elasticity_matrix(3,4) = 0.0; elasticity_matrix(3,5) = 0.0;
  elasticity_matrix(4,0) = 0.0; elasticity_matrix(4,1) = 0.0; elasticity_matrix(4,2) = 0.0; elasticity_matrix(4,3) = 0.0; elasticity_matrix(4,4) = c44; elasticity_matrix(4,5) = 0.0;
  elasticity_matrix(5,0) = 0.0; elasticity_matrix(5,1) = 0.0; elasticity_matrix(5,2) = 0.0; elasticity_matrix(5,3) = 0.0; elasticity_matrix(5,4) = 0.0; elasticity_matrix(5,5) = c44;

  double gphi = (1.0-phi)*(1.0-phi) + 1.0e-6;

  elasticity_matrix.Scale(gphi);

  /*
  char jobz = 'V';
  char uplo = 'U';

  double * ap = new double[6];

  ap[0] = epsilon(0); ap[1] = epsilon(5); ap[3] = epsilon(4);
                      ap[2] = epsilon(1); ap[4] = epsilon(3);
                                          ap[5] = epsilon(2);

  double * WORK = new double[6];
  int INFO;

  Epetra_SerialDenseVector eigv(3);
  Epetra_SerialDenseVector eigw(9);

  Lapack->SPEV(jobz, uplo, 3, ap, &eigv[0], &eigw[0], 3, WORK, &INFO);

  Epetra_SerialDenseVector n1(3), n2(3), n3(3);
  n1(0) = eigw(0); n1(1) = eigw(1); n1(2) = eigw(2);
  n2(0) = eigw(3); n2(1) = eigw(4); n2(2) = eigw(5);
  n3(0) = eigw(6); n3(1) = eigw(7); n3(2) = eigw(8);


  Epetra_SerialDenseVector eye(6);
  Epetra_SerialDenseVector m1(6), m2(6), m3(6);
  Epetra_SerialDenseMatrix Ep(6,6), En(6,6);

  Epetra_SerialDenseMatrix elasticity_neg(6,6);

  Epetra_SerialDenseMatrix G12(6,6), G13(6,6), G23(6,6), G(6,6);

  m1(0) = n1(0)*n1(0);
  m1(1) = n1(1)*n1(1);
  m1(2) = n1(2)*n1(2);
  m1(3) = n1(1)*n1(2);
  m1(4) = n1(0)*n1(2);
  m1(5) = n1(0)*n1(1);

  m2(0) = n2(0)*n2(0);
  m2(1) = n2(1)*n2(1);
  m2(2) = n2(2)*n2(2);
  m2(3) = n2(1)*n2(2);
  m2(4) = n2(0)*n2(2);
  m2(5) = n2(0)*n2(1);

  m3(0) = n3(0)*n3(0);
  m3(1) = n3(1)*n3(1);
  m3(2) = n3(2)*n3(2);
  m3(3) = n3(1)*n3(2);
  m3(4) = n3(0)*n3(2);
  m3(5) = n3(0)*n3(1);

  G12(0,0) = 4.0*m1(0)*m2(0); G12(0,1) = 4.0*m1(5)*m2(5); G12(0,2) = 4.0*m1(4)*m2(4); G12(0,3) = 2.0*m1(4)*m2(5) + 2.0*m1(5)*m2(4); G12(0,4) = 2.0*m1(0)*m2(4) + 2.0*m1(4)*m2(0); G12(0,5) = 2.0*m1(0)*m2(5) + 2.0*m1(5)*m2(0);
  G12(1,0) = 4.0*m1(5)*m2(5); G12(1,1) = 4.0*m1(1)*m2(1); G12(1,2) = 4.0*m1(3)*m2(3); G12(1,3) = 2.0*m1(1)*m2(3) + 2.0*m1(3)*m2(1); G12(1,4) = 2.0*m1(3)*m2(5) + 2.0*m1(5)*m2(3); G12(1,5) = 2.0*m1(1)*m2(5) + 2.0*m1(5)*m2(1);
  G12(2,0) = 4.0*m1(4)*m2(4); G12(2,1) = 4.0*m1(3)*m2(3); G12(2,2) = 4.0*m1(2)*m2(2); G12(2,3) = 2.0*m1(2)*m2(3) + 2.0*m1(3)*m2(2); G12(2,4) = 2.0*m1(2)*m2(4) + 2.0*m1(4)*m2(2); G12(2,5) = 2.0*m1(3)*m2(4) + 2.0*m1(4)*m2(3);
  G12(3,0) = 2.0*m1(4)*m2(5) + 2.0*m1(5)*m2(4); G12(3,1) = 2.0*m1(1)*m2(3) + 2.0*m1(3)*m2(1); G12(3,2) = 2.0*m1(2)*m2(3) + 2.0*m1(3)*m2(2); G12(3,3) = m1(1)*m2(2) + m1(2)*m2(1) + 2.0*m1(3)*m2(3); G12(3,4) = m1(2)*m2(5) + m1(3)*m2(4) + m1(4)*m2(3) + m1(5)*m2(2); G12(3,5) = m1(1)*m2(4) + m1(4)*m2(1) + m1(3)*m2(5) + m1(5)*m2(3);
  G12(4,0) = 2.0*m1(0)*m2(4) + 2.0*m1(4)*m2(0); G12(4,1) = 2.0*m1(3)*m2(5) + 2.0*m1(5)*m2(3); G12(4,2) = 2.0*m1(2)*m2(4) + 2.0*m1(4)*m2(2); G12(4,3) = m1(2)*m2(5) + m1(3)*m2(4) + m1(4)*m2(3) + m1(5)*m2(2); G12(4,4) = m1(0)*m2(2) + m1(2)*m2(0) + 2.0*m1(4)*m2(4); G12(4,5) = m1(0)*m2(3) + m1(3)*m2(0) + m1(4)*m2(5) + m1(5)*m2(4);
  G12(5,0) = 2.0*m1(0)*m2(5) + 2.0*m1(5)*m2(0); G12(5,1) = 2.0*m1(1)*m2(5) + 2.0*m1(5)*m2(1); G12(5,2) = 2.0*m1(3)*m2(4) + 2.0*m1(4)*m2(3); G12(5,3) = m1(1)*m2(4) + m1(4)*m2(1) + m1(3)*m2(5) + m1(5)*m2(3); G12(5,4) = m1(0)*m2(3) + m1(3)*m2(0) + m1(4)*m2(5) + m1(5)*m2(4); G12(5,5) = m1(0)*m2(1) + m1(1)*m2(0) + 2.0*m1(5)*m2(5);

  G13(0,0) = 4.0*m1(0)*m3(0); G13(0,1) = 4.0*m1(5)*m3(5); G13(0,2) = 4.0*m1(4)*m3(4); G13(0,3) = 2.0*m1(4)*m3(5) + 2.0*m1(5)*m3(4); G13(0,4) = 2.0*m1(0)*m3(4) + 2.0*m1(4)*m3(0); G13(0,5) = 2.0*m1(0)*m3(5) + 2.0*m1(5)*m3(0);
  G13(1,0) = 4.0*m1(5)*m3(5); G13(1,1) = 4.0*m1(1)*m3(1); G13(1,2) = 4.0*m1(3)*m3(3); G13(1,3) = 2.0*m1(1)*m3(3) + 2.0*m1(3)*m3(1); G13(1,4) = 2.0*m1(3)*m3(5) + 2.0*m1(5)*m3(3); G13(1,5) = 2.0*m1(1)*m3(5) + 2.0*m1(5)*m3(1);
  G13(2,0) = 4.0*m1(4)*m3(4); G13(2,1) = 4.0*m1(3)*m3(3); G13(2,2) = 4.0*m1(2)*m3(2); G13(2,3) = 2.0*m1(2)*m3(3) + 2.0*m1(3)*m3(2); G13(2,4) = 2.0*m1(2)*m3(4) + 2.0*m1(4)*m3(2); G13(2,5) = 2.0*m1(3)*m3(4) + 2.0*m1(4)*m3(3);
  G13(3,0) = 2.0*m1(4)*m3(5) + 2.0*m1(5)*m3(4); G13(3,1) = 2.0*m1(1)*m3(3) + 2.0*m1(3)*m3(1); G13(3,2) = 2.0*m1(2)*m3(3) + 2.0*m1(3)*m3(2); G13(3,3) = m1(1)*m3(2) + m1(2)*m3(1) + 2.0*m1(3)*m3(3); G13(3,4) = m1(2)*m3(5) + m1(3)*m3(4) + m1(4)*m3(3) + m1(5)*m3(2); G13(3,5) = m1(1)*m3(4) + m1(4)*m3(1) + m1(3)*m3(5) + m1(5)*m3(3);
  G13(4,0) = 2.0*m1(0)*m3(4) + 2.0*m1(4)*m3(0); G13(4,1) = 2.0*m1(3)*m3(5) + 2.0*m1(5)*m3(3); G13(4,2) = 2.0*m1(2)*m3(4) + 2.0*m1(4)*m3(2); G13(4,3) = m1(2)*m3(5) + m1(3)*m3(4) + m1(4)*m3(3) + m1(5)*m3(2); G13(4,4) = m1(0)*m3(2) + m1(2)*m3(0) + 2.0*m1(4)*m3(4); G13(4,5) = m1(0)*m3(3) + m1(3)*m3(0) + m1(4)*m3(5) + m1(5)*m3(4);
  G13(5,0) = 2.0*m1(0)*m3(5) + 2.0*m1(5)*m3(0); G13(5,1) = 2.0*m1(1)*m3(5) + 2.0*m1(5)*m3(1); G13(5,2) = 2.0*m1(3)*m3(4) + 2.0*m1(4)*m3(3); G13(5,3) = m1(1)*m3(4) + m1(4)*m3(1) + m1(3)*m3(5) + m1(5)*m3(3); G13(5,4) = m1(0)*m3(3) + m1(3)*m3(0) + m1(4)*m3(5) + m1(5)*m3(4); G13(5,5) = m1(0)*m3(1) + m1(1)*m3(0) + 2.0*m1(5)*m3(5);

  G23(0,0) = 4.0*m2(0)*m3(0); G23(0,1) = 4.0*m2(5)*m3(5); G23(0,2) = 4.0*m2(4)*m3(4); G23(0,3) = 2.0*m2(4)*m3(5) + 2.0*m2(5)*m3(4); G23(0,4) = 2.0*m2(0)*m3(4) + 2.0*m2(4)*m3(0); G23(0,5) = 2.0*m2(0)*m3(5) + 2.0*m2(5)*m3(0);
  G23(1,0) = 4.0*m2(5)*m3(5); G23(1,1) = 4.0*m2(1)*m3(1); G23(1,2) = 4.0*m2(3)*m3(3); G23(1,3) = 2.0*m2(1)*m3(3) + 2.0*m2(3)*m3(1); G23(1,4) = 2.0*m2(3)*m3(5) + 2.0*m2(5)*m3(3); G23(1,5) = 2.0*m2(1)*m3(5) + 2.0*m2(5)*m3(1);
  G23(2,0) = 4.0*m2(4)*m3(4); G23(2,1) = 4.0*m2(3)*m3(3); G23(2,2) = 4.0*m2(2)*m3(2); G23(2,3) = 2.0*m2(2)*m3(3) + 2.0*m2(3)*m3(2); G23(2,4) = 2.0*m2(2)*m3(4) + 2.0*m2(4)*m3(2); G23(2,5) = 2.0*m2(3)*m3(4) + 2.0*m2(4)*m3(3);
  G23(3,0) = 2.0*m2(4)*m3(5) + 2.0*m2(5)*m3(4); G23(3,1) = 2.0*m2(1)*m3(3) + 2.0*m2(3)*m3(1); G23(3,2) = 2.0*m2(2)*m3(3) + 2.0*m2(3)*m3(2); G23(3,3) = m2(1)*m3(2) + m2(2)*m3(1) + 2.0*m2(3)*m3(3); G23(3,4) = m2(2)*m3(5) + m2(3)*m3(4) + m2(4)*m3(3) + m2(5)*m3(2); G23(3,5) = m2(1)*m3(4) + m2(4)*m3(1) + m2(3)*m3(5) + m2(5)*m3(3);
  G23(4,0) = 2.0*m2(0)*m3(4) + 2.0*m2(4)*m3(0); G23(4,1) = 2.0*m2(3)*m3(5) + 2.0*m2(5)*m3(3); G23(4,2) = 2.0*m2(2)*m3(4) + 2.0*m2(4)*m3(2); G23(4,3) = m2(2)*m3(5) + m2(3)*m3(4) + m2(4)*m3(3) + m2(5)*m3(2); G23(4,4) = m2(0)*m3(2) + m2(2)*m3(0) + 2.0*m2(4)*m3(4); G23(4,5) = m2(0)*m3(3) + m2(3)*m3(0) + m2(4)*m3(5) + m2(5)*m3(4);
  G23(5,0) = 2.0*m2(0)*m3(5) + 2.0*m2(5)*m3(0); G23(5,1) = 2.0*m2(1)*m3(5) + 2.0*m2(5)*m3(1); G23(5,2) = 2.0*m2(3)*m3(4) + 2.0*m2(4)*m3(3); G23(5,3) = m2(1)*m3(4) + m2(4)*m3(1) + m2(3)*m3(5) + m2(5)*m3(3); G23(5,4) = m2(0)*m3(3) + m2(3)*m3(0) + m2(4)*m3(5) + m2(5)*m3(4); G23(5,5) = m2(0)*m3(1) + m2(1)*m3(0) + 2.0*m2(5)*m3(5);

  double tol = 1E-10;
  if (std::fabs(eigv[0]-eigv[1])>tol){
    G = G12;
    G.Scale( (std::max(0.0,eigv[0]) - std::max(0.0,eigv[1]))/(2.0*(eigv[0]-eigv[1])) );
    Ep = G;

    G = G12;
    G.Scale( (std::min(0.0,eigv[0]) - std::min(0.0,eigv[1]))/(2.0*(eigv[0]-eigv[1])) );
    En = G;
  }
  else{
    G = G12;
    G.Scale(0.5);
    Ep = G;
    En = G;
  }
  if (std::fabs(eigv[0]-eigv[2])>tol){
    G = G13;
    G.Scale( (std::max(0.0,eigv[0]) - std::max(0.0,eigv[2]))/(2.0*(eigv[0]-eigv[2])) );
    Ep += G;

    G = G13;
    G.Scale( (std::min(0.0,eigv[0]) - std::min(0.0,eigv[2]))/(2.0*(eigv[0]-eigv[2])) );
    En += G;
  }
  else{
    G = G13;
    G.Scale(0.5);
    Ep += G;
    En += G;
  }
  if (std::fabs(eigv[1]-eigv[2])>tol){
    G = G23;
    G.Scale( (std::max(0.0,eigv[1]) - std::max(0.0,eigv[2]))/(2.0*(eigv[1]-eigv[2])) );
    Ep += G;

    G = G23;
    G.Scale( (std::min(0.0,eigv[1]) - std::min(0.0,eigv[2]))/(2.0*(eigv[1]-eigv[2])) );
    En += G;
  }
  else{
    G = G23;
    G.Scale(0.5);
    Ep += G;
    En += G;
  }

  if (eigv(0)>0.0){
    Ep.Multiply('N','T',1.0,m1,m1,1.0);
  }
  if (eigv(0)<0.0){
    En.Multiply('N','T',1.0,m1,m1,1.0);
  }

  if (eigv(1)>0.0){
    Ep.Multiply('N','T',1.0,m2,m2,1.0);
  }
  if (eigv(1)<0.0){
    En.Multiply('N','T',1.0,m2,m2,1.0);
  }

  if (eigv(2)>0.0){
    Ep.Multiply('N','T',1.0,m3,m3,1.0);
  }
  if (eigv(2)<0.0){
    En.Multiply('N','T',1.0,m3,m3,1.0);
  }

  double tr_e = epsilon(0) + epsilon(1) + epsilon(2);
  eye(0) = 1.0; eye(1) = 1.0; eye(2) = 1.0; eye(3) = 0.0; eye(4) = 0.0; eye(5) = 0.0;

  if (tr_e > 0.0){
    elasticity_matrix.Multiply('N','T',lambda,eye,eye,0.0);
  }
  if (tr_e < 0.0){
    elasticity_neg.Multiply('N','T',lambda,eye,eye,0.0);
  }

  Ep.Scale(2.0*mu);
  En.Scale(2.0*mu);

  elasticity_matrix += Ep;
  elasticity_neg    += En;

  elasticity_matrix.Scale(gphi);

  elasticity_matrix += elasticity_neg;
  */
}

void TPF::elasticity::setup_dirichlet_conditions(){

  // example of how the setup should be implemented
  // this should be simplified or hidden in next versions

  n_bc_dof = 0;
  double z;
  int node;

  // get the total number of prescribed dof
  for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
    node = Mesh->local_nodes[i];
    z    = Mesh->nodes_coord[3*node+2]; // where do you want
    if(z<=1.0e-6 && z>=-1.0e-6){        // to apply Dirichlet conditions ?
        n_bc_dof+=3;                    // how many dof at this location ?
    }
    if(z<=10+1.0e-6 && z>=10-1.0e-6){   // how many dof at this location ?
        n_bc_dof+=1;
    }
  }

  // store the prescibred dof into a vector dof_on_boundary
  int indbc = 0;
  dof_on_boundary = new int [n_bc_dof];
  for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
      node = Mesh->local_nodes[inode];
      z    = Mesh->nodes_coord[3*node+2];
      if (z<=1.0e-6 && z>=-1.0e-6){
          dof_on_boundary[indbc+0] = 3*inode+0;
          dof_on_boundary[indbc+1] = 3*inode+1;
          dof_on_boundary[indbc+2] = 3*inode+2;
          indbc+=3;
      }
      if (z<=10+1.0e-6 && z>=10-1.0e-6){
          dof_on_boundary[indbc+0] = 3*inode+2;
          indbc+=1;
      }
  }

}

void TPF::elasticity::apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){

  // example of how dirichlet conditions should be applied
  // this should be simplified or hidden in next versions

  Epetra_MultiVector v(*Mesh->StandardMapU,true);
  v.PutScalar(0.0);

  int node;
  double z;
  for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
      node = Mesh->local_nodes[inode];
      z    = Mesh->nodes_coord[3*node+2];
      if (z<=10+1.0e-6 && z>=10-1.0e-6){
          v[0][Mesh->StandardMapU->LID(3*node+2)] = displacement;
      }
  }

  Epetra_MultiVector rhs_dir(*Mesh->StandardMapU,true);
  K.Apply(v,rhs_dir);
  F.Update(-1.0,rhs_dir,1.0);

  for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
      node = Mesh->local_nodes[inode];
      z    = Mesh->nodes_coord[3*node+2];
      if (z<=1.0e-6 && z>=-1.0e-6){
          F[0][Mesh->StandardMapU->LID(3*node+0)] = 0.0;
          F[0][Mesh->StandardMapU->LID(3*node+1)] = 0.0;
          F[0][Mesh->StandardMapU->LID(3*node+2)] = 0.0;
      }
      if (z<=10+1.0e-6 && z>=10-1.0e-6){
          F[0][Mesh->StandardMapU->LID(3*node+2)] = displacement;
      }
  }
  //}
  ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
}

void TPF::elasticity::updateDamageHistory(Epetra_Vector & damageHistory, Epetra_Vector & u){

  Epetra_Vector v(*Mesh->OverlapMapU);
  v.Import(u, *Mesh->ImportToOverlapMapU, Insert);

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
      dx_shape_functions(inode,0) = Mesh->DX_N_cells(inode,elid);
      dx_shape_functions(inode,1) = Mesh->DY_N_cells(inode,elid);
      dx_shape_functions(inode,2) = Mesh->DZ_N_cells(inode,elid);
      for (unsigned int ddl=0; ddl<3; ++ddl){
        id = 3*node+ddl;
        cells_u(3*inode+ddl) = v[Mesh->OverlapMapU->LID(id)];
      }
    }

    compute_B_matrices(dx_shape_functions,matrix_B);
    epsilon.Multiply('N', 'N', 1.0, matrix_B,cells_u, 0.0);

    trepsilon  = epsilon(0) + epsilon(1) + epsilon(2);
    trepsilon2 = epsilon(0)*epsilon(0) + epsilon(1)*epsilon(1) + epsilon(2)*epsilon(2) +
                 0.5*epsilon(3)*epsilon(3) + 0.5*epsilon(4)*epsilon(4) + 0.5*epsilon(5)*epsilon(5);

    history   = damageHistory[elid];
    potential = (lambda/2.0)*trepsilon*trepsilon + mu*trepsilon2;

    if (potential>history){
      damageHistory[elid] = potential;
    }
  }

}
