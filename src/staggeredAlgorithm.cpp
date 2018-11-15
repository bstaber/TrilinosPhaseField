#include "staggeredAlgorithm.hpp"

staggeredAlgorithm::staggeredAlgorithm(Epetra_Comm & comm, mesh & mesh_, damageBVP & damageInterface_, elasticBVP & elasticInterface_){

  Comm = &comm;
  Mesh = &mesh_;

  damageInterface = &damageInterface_;
  elasticInterface = &elasticInterface_;
}

staggeredAlgorithm::~staggeredAlgorithm(){
}

Epetra_Map staggeredAlgorithm::constructGaussMap(){
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

void staggeredAlgorithm::staggeredAlgorithmDirichletBC(Teuchos::ParameterList & ParametersList, bool print){

  double delta_u  = Teuchos::getParameter<double>(ParametersList.sublist("Elasticity"), "delta_u");
  int n_steps     = Teuchos::getParameter<int>(ParametersList.sublist("Elasticity"), "n_steps");

  Epetra_Time Time(*Comm);

  Epetra_FECrsMatrix matrix_d(Copy,*damageInterface->FEGraph);
  Epetra_FECrsMatrix matrix_u(Copy,*elasticInterface->FEGraph);

  Epetra_FEVector    rhs_d(*damageInterface->StandardMap);
  Epetra_FEVector    rhs_u(*elasticInterface->StandardMap);

  Epetra_Vector      lhs_d(*damageInterface->StandardMap);
  Epetra_Vector      lhs_u(*elasticInterface->StandardMap);

  Epetra_Vector      phi(*damageInterface->OverlapMap);

  Epetra_Map GaussMap = constructGaussMap();
  Epetra_Vector damageHistory(GaussMap);

  if (Comm->MyPID()==0){
    std::cout << "step" << std::setw(15) << "cpu_time (s)" << "\n";
  }

  double bc_disp = 0.0;
  damageHistory.PutScalar(0.0);
  for (int n=0; n<n_steps; ++n){

    Time.ResetStartTime();
    damageInterface->solve(ParametersList.sublist("Damage"), matrix_d, lhs_d, rhs_d, damageHistory, GaussMap);

    phi.Import(lhs_d, *damageInterface->ImportToOverlapMap, Insert);

    bc_disp = (double(n)+1.0)*delta_u;
    elasticInterface->computeDisplacement(ParametersList.sublist("Elasticity"), matrix_u, lhs_u, rhs_u,
                                                                                phi, *damageInterface->OverlapMap,
                                                                                bc_disp);
    elasticInterface->updateDamageHistory(damageHistory, lhs_u, GaussMap);

    if (Comm->MyPID()==0){
      std::cout << n << std::setw(15) << Time.ElapsedTime() << "\n";
    }
  }

  if (print){
    std::string printpath = Teuchos::getParameter<std::string>(ParametersList.sublist("Mesh"),"printpath");
    std::string dispfile  = printpath + "displacement" + std::to_string(int(n_steps)) + ".mtx";
    std::string damgfile  = printpath + "damage"       + std::to_string(int(n_steps)) + ".mtx";
    int error_u = elasticInterface->print_solution(lhs_u, dispfile);
    int error_d = damageInterface->print_solution(lhs_d, damgfile);
  }

}
