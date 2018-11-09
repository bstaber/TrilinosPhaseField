#include "staggeredAlgorithm.hpp"

TPF::staggeredAlgorithm::staggeredAlgorithm(Epetra_Comm & comm, mesh & mesh_)
: elasticity(comm, mesh_){

  Lapack                        = new Epetra_LAPACK();
  displacementSolutionOverlaped = new Epetra_Vector(*Mesh->OverlapMapU);
  damageSolutionOverlaped       = new Epetra_Vector(*Mesh->OverlapMapD);
}

TPF::staggeredAlgorithm::~staggeredAlgorithm(){

}

void TPF::staggeredAlgorithm::staggeredAlgorithmDirichletBC(Teuchos::ParameterList & ParametersList, bool print){

  gc = Teuchos::getParameter<double>(ParametersList.sublist("Damage"), "gc");
  lc = Teuchos::getParameter<double>(ParametersList.sublist("Damage"), "lc");

  Teuchos::RCP<damage> phaseFieldBVP = Teuchos::rcp(new damage(*Comm, *Mesh, gc, lc));

  double delta_u  = Teuchos::getParameter<double>(ParametersList.sublist("Elasticity"), "delta_u");
  int n_steps     = Teuchos::getParameter<int>(ParametersList.sublist("Elasticity"), "n_steps");

  Epetra_Time Time(*Comm);

  Epetra_FECrsMatrix matrix_d(Copy,*Mesh->FEGraphD);
  Epetra_FECrsMatrix matrix_u(Copy,*Mesh->FEGraphU);

  Epetra_FEVector    rhs_d(*Mesh->StandardMapD);
  Epetra_FEVector    rhs_u(*Mesh->StandardMapU);

  Epetra_Vector      lhs_d(*Mesh->StandardMapD);
  Epetra_Vector      lhs_u(*Mesh->StandardMapU);

  Epetra_Map ElementMap(-1, Mesh->n_local_cells, &Mesh->local_cells[0], 0, *Comm);
  Epetra_Vector damageHistory(ElementMap);

  if (Comm->MyPID()==0){
    std::cout << "step" << std::setw(15) << "cpu_time (s)" << "\n";
  }

  double bc_disp = 0.0;
  damageHistory.PutScalar(0.0);
  damageSolutionOverlaped->PutScalar(0.0);
  displacementSolutionOverlaped->PutScalar(0.0);

  for (int n=0; n<n_steps; ++n){

    Time.ResetStartTime();

    bc_disp = (double(n)+1.0)*delta_u;

    solve_u(matrix_u, rhs_u, lhs_u, ParametersList, bc_disp);

    updateDamageHistory(damageHistory);

    phaseFieldBVP->solve_d(matrix_d, rhs_d, lhs_d, ParametersList, damageHistory);

    damageSolutionOverlaped->Import(lhs_d, *Mesh->ImportToOverlapMapD, Insert);
    displacementSolutionOverlaped->Import(lhs_u, *Mesh->ImportToOverlapMapU, Insert);

    if (Comm->MyPID()==0){
      std::cout << n << std::setw(15) << Time.ElapsedTime() << "\n";
    }

    if (print){
      std::string dispfile = "/home/s/staber/Trilinos_results/examples/phasefield/displacement" + std::to_string(int(n)) + ".mtx";
      std::string damgfile = "/home/s/staber/Trilinos_results/examples/phasefield/damage"       + std::to_string(int(n)) + ".mtx";
      int error_u = print_solution(lhs_u, dispfile);
      int error_d = phaseFieldBVP->print_solution(lhs_d, damgfile);
    }
    
  }

}

void TPF::staggeredAlgorithm::updateDamageHistory(Epetra_Vector & damageHistory){

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
        cells_u(3*inode+ddl) = displacementSolutionOverlaped[0][Mesh->OverlapMapU->LID(id)];
      }
    }

    compute_B_matrices(dx_shape_functions,matrix_B);
    epsilon.Multiply('N','N',1.0,matrix_B,cells_u,0.0);

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

void TPF::staggeredAlgorithm::get_elasticity_tensor(Epetra_SerialDenseMatrix & elasticity_matrix, Epetra_SerialDenseVector & epsilon, double & phi){

  double gphi = (1.0-phi)*(1.0-phi) + 1.0e-6;

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

  Epetra_SerialDenseMatrix elasticity_pos(6,6);
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
    elasticity_pos.Multiply('N','T',lambda,eye,eye,0.0);
  }
  if (tr_e < 0.0){
    elasticity_neg.Multiply('N','T',lambda,eye,eye,0.0);
  }

  Ep.Scale(2.0*mu);
  En.Scale(2.0*mu);

  elasticity_pos += Ep;
  elasticity_neg += En;

  elasticity_pos.Scale(gphi);

}
