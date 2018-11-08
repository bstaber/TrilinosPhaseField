#include "staggeredAlgorithm.hpp"

TPF::staggeredAlgorithm::staggeredAlgorithm(Epetra_Comm & comm, mesh & mesh_) : elasticity(comm, mesh_){
  Mesh = &mesh_;
  Comm = &comm;

  phaseFieldBVP = Teuchos::rcp(new damage(comm, *Mesh, gc, lc));
}

TPF::staggeredAlgorithm::~staggeredAlgorithm(){
}

void TPF::staggeredAlgorithm::staggeredAlgorithmDirichletBC(Teuchos::ParameterList & ParametersList, bool print){

  double delta_u  = Teuchos::getParameter<double>(ParametersList.sublist("Elasticity"), "delta_u");
  int n_steps     = Teuchos::getParameter<int>(ParametersList.sublist("Elasticity"), "n_steps");

  Epetra_Time Time(*Comm);

  Epetra_FECrsMatrix matrix_d(Copy,*Mesh->FEGraphD);
  Epetra_FECrsMatrix matrix_u(Copy,*Mesh->FEGraphU);

  Epetra_FEVector    rhs_d(*Mesh->StandardMapD);
  Epetra_FEVector    rhs_u(*Mesh->StandardMapU);

  Epetra_Vector      lhs_d(*Mesh->StandardMapD);
  Epetra_Vector      lhs_u(*Mesh->StandardMapU);

  /*
  Epetra_Map GaussMap = constructGaussMap();
  Epetra_Vector damageHistory(GaussMap);
  */

  if (Comm->MyPID()==0){
    std::cout << "step" << std::setw(15) << "cpu_time (s)" << "\n";
  }

  double bc_disp = 0.0;
  //damageHistory.PutScalar(0.0);
  for (int n=0; n<n_steps; ++n){

    Time.ResetStartTime();

    bc_disp = (double(n)+1.0)*delta_u;

    //damageSolutionOverlaped->Import(lhs_d, *damageInterface->ImportToOverlapMap, Insert);

    //updateDamageHistory(damageHistory, lhs_u, GaussMap);

    //phaseFieldBVP->solve(ParametersList.sublist("Damage"), matrix_d, lhs_d, rhs_d, damageHistory);

    //computeDisplacement(ParametersList.sublist("Elasticity"), matrix_u, lhs_u, rhs_u, bc_disp);

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

/*
void TPF::staggeredAlgorithm::computeDisplacement(Teuchos::ParameterList & Parameters,
                                                    Epetra_FECrsMatrix & matrix, Epetra_Vector & lhs, Epetra_FEVector & rhs,
                                                    double & bc_disp){

  rhs.PutScalar(0.0);
  lhs.PutScalar(0.0);

  assemblePureDirichlet_homogeneousForcing(matrix);
  apply_dirichlet_conditions(matrix, rhs, bc_disp);

  int max_iter = Teuchos::getParameter<int>(Parameters.sublist("Aztec"), "AZ_max_iter");
  double tol   = Teuchos::getParameter<double>(Parameters.sublist("Aztec"), "AZ_tol");

  Epetra_LinearProblem problem(&matrix, &lhs, &rhs);

  AztecOO solver(problem);
  solver.SetParameters(Parameters.sublist("Aztec"));
  solver.Iterate(max_iter, tol);
}

void TPF::staggeredAlgorithm::updateDamageHistory(Epetra_Vector & damageHistory,
                                                    Epetra_Vector & displacement){

  Epetra_Vector u(*OverlapMap);
  u.Import(displacement, *ImportToOverlapMap, Insert);

  //Mesh->update_store_feinterp_cells(u, *OverlapMap);

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

Epetra_Map TPF::staggeredAlgorithm::constructGaussMap(){
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

void TPF::staggeredAlgorithm::get_elasticity_tensor(unsigned int & e_lid,
                                                           unsigned int & gp,
                                                           Epetra_SerialDenseMatrix & tangent_matrix){
  int n_gauss_points = Mesh->n_gauss_cells;
  int e_gid = Mesh->local_cells[e_lid];
  int node;
  double xi = Mesh->xi_cells[gp]; double eta = Mesh->eta_cells[gp]; double zeta = Mesh->zeta_cells[gp];

  Epetra_SerialDenseVector shape_functions(Mesh->el_type);
  switch (Mesh->el_type){
      case 4:
          tetra4::shape_functions(shape_functions, xi, eta, zeta);
          break;
      case 8:
          hexa8::shape_functions(shape_functions, xi, eta, zeta);
          break;
      case 10:
          tetra10::shape_functions(shape_functions, xi, eta, zeta);
          break;
  }

  double d = 0.0;
  for (unsigned int j=0; j<Mesh->el_type; ++j){
      node = Mesh->cells_nodes[Mesh->el_type*e_gid+j];
      d += shape_functions(j) * damageSolutionOverlaped[0][damageInterface->OverlapMap->LID(node)];
  }

  double g = (1.0-d)*(1.0-d) + 1.0e-6;
  tangent_matrix = elasticity;
  tangent_matrix.Scale(g);
}

void TPF::staggeredAlgorithm::get_elasticity_tensor_for_recovery(unsigned int & e_lid,
                                                                        Epetra_SerialDenseMatrix & tangent_matrix){
  std::cout << "Not using this method yet.\n";
}
*/
