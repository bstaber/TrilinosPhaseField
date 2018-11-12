#include "staggeredAlgorithm.hpp"

TPF::staggeredAlgorithm::staggeredAlgorithm(Epetra_Comm & comm, mesh & mesh_)
: elasticity(comm, mesh_){

}

TPF::staggeredAlgorithm::~staggeredAlgorithm(){

}

void TPF::staggeredAlgorithm::staggeredAlgorithmDirichletBC(Teuchos::ParameterList & ParametersList, bool print){

  std::string path = Teuchos::getParameter<std::string>(ParametersList.sublist("Mesh"), "path_to_results");

  gc = Teuchos::getParameter<double>(ParametersList.sublist("Damage"), "gc");
  lc = Teuchos::getParameter<double>(ParametersList.sublist("Damage"), "lc");
  Teuchos::RCP<damage> phaseFieldBVP = Teuchos::rcp(new damage(*Comm, *Mesh, gc, lc));

  double E  = Teuchos::getParameter<double>(ParametersList.sublist("Elasticity"), "young");
  double nu = Teuchos::getParameter<double>(ParametersList.sublist("Elasticity"), "poisson");

  lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
  mu     = E/(2.0*(1.0+nu));

  double delta_u  = Teuchos::getParameter<double>(ParametersList.sublist("Elasticity"), "delta_u");
  int n_steps     = Teuchos::getParameter<int>(ParametersList.sublist("Elasticity"), "n_steps");

  Epetra_Time Time(*Comm);

  Epetra_FECrsMatrix matrix_d(Copy,*Mesh->FEGraphD);
  Epetra_FECrsMatrix matrix_u(Copy,*Mesh->FEGraphU);

  Epetra_FEVector    rhs_d(*Mesh->StandardMapD);
  Epetra_FEVector    rhs_u(*Mesh->StandardMapU);

  Epetra_Vector      lhs_d(*Mesh->StandardMapD);
  Epetra_Vector      lhs_u(*Mesh->StandardMapU);

  Epetra_Vector phi(*Mesh->OverlapMapD);

  Epetra_Map ElementMap(-1, Mesh->n_local_cells, &Mesh->local_cells[0], 0, *Comm);
  Epetra_Vector damageHistory(ElementMap);

  double bc_disp = 0.0;
  damageHistory.PutScalar(0.0);

  lhs_d.PutScalar(0.0);
  lhs_u.PutScalar(0.0);

  if (Comm->MyPID()==0){
    std::cout << "step" << std::setw(15) << "cpu_time (s)" << "\n";
  }

  Epetra_Vector u(*Mesh->StandardMapU);
  for (int n=0; n<n_steps; ++n){

    Time.ResetStartTime();

    bc_disp = (double(n)+1.0)*delta_u;

    phi.Import(lhs_d, *Mesh->ImportToOverlapMapD, Insert);
    u = lhs_u;

    lhs_u = solve_u(matrix_u, rhs_u, u, ParametersList, bc_disp, phi);


    updateDamageHistory(damageHistory, lhs_u);

    phaseFieldBVP->solve_d(matrix_d, rhs_d, lhs_d, ParametersList, damageHistory);

    if (Comm->MyPID()==0){
      std::cout << n << std::setw(15) << Time.ElapsedTime() << "\n";
    }

    /*
    std::string dispfile = path + "displacement" + std::to_string(int(n)) + ".mtx";
    std::string damgfile = path + "damage"       + std::to_string(int(n)) + ".mtx";
    int error_u = print_solution(lhs_u, dispfile);
    int error_d = phaseFieldBVP->print_solution(lhs_d, damgfile);
    */
  }

}

void TPF::staggeredAlgorithm::updateDamageHistory(Epetra_Vector & damageHistory, Epetra_Vector & u){

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
