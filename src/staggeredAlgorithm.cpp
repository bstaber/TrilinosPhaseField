#include "staggeredAlgorithm.hpp"

TPF::staggeredAlgorithm::staggeredAlgorithm(Epetra_Comm & comm, mesh & mesh_){
  Mesh = &mesh_;
  Comm = &comm;
}

TPF::staggeredAlgorithm::~staggeredAlgorithm(){

}

void TPF::staggeredAlgorithm::staggeredAlgorithmDirichletBC(Teuchos::ParameterList & ParametersList, bool print){

  std::string path = Teuchos::getParameter<std::string>(ParametersList.sublist("Mesh"), "path_to_results");

  double gc = Teuchos::getParameter<double>(ParametersList.sublist("Damage"), "gc");
  double lc = Teuchos::getParameter<double>(ParametersList.sublist("Damage"), "lc");

  double E  = Teuchos::getParameter<double>(ParametersList.sublist("Elasticity"), "young");
  double nu = Teuchos::getParameter<double>(ParametersList.sublist("Elasticity"), "poisson");
  double lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
  double mu     = E/(2.0*(1.0+nu));

  Teuchos::RCP<damage> phaseFieldBVP = Teuchos::rcp(new damage(*Comm, *Mesh, gc, lc));
  Teuchos::RCP<elasticity> elasticityBVP = Teuchos::rcp(new elasticity(*Comm, *Mesh, lambda, mu));

  double delta_u = Teuchos::getParameter<double>(ParametersList.sublist("Elasticity"), "delta_u");
  int n_steps    = Teuchos::getParameter<int>(ParametersList.sublist("Elasticity"), "n_steps");

  Epetra_Time Time(*Comm);

  Epetra_FECrsMatrix matrix_d(Copy,*phaseFieldBVP->FEGraph);
  Epetra_FECrsMatrix matrix_u(Copy,*elasticityBVP->FEGraph);

  Epetra_FEVector    rhs_d(*phaseFieldBVP->StandardMap);
  Epetra_FEVector    rhs_u(*elasticityBVP->StandardMap);

  Epetra_Vector      lhs_d(*phaseFieldBVP->StandardMap);
  Epetra_Vector      lhs_u(*elasticityBVP->StandardMap);

  Epetra_Vector      w(*phaseFieldBVP->OverlapMap);
  Epetra_Vector      v(*elasticityBVP->OverlapMap);


  Epetra_Map ElementMap(-1, Mesh->n_local_cells, &Mesh->local_cells[0], 0, *Comm);
  Epetra_Vector damageHistory(ElementMap);

  double bc_disp = 0.0;

  damageHistory.PutScalar(0.0);
  lhs_d.PutScalar(0.0);
  lhs_u.PutScalar(0.0);

  if (Comm->MyPID()==0){
    std::cout << "step" << std::setw(15) << "cpu_time (s)" << "\n";
  }

  for (int n=0; n<n_steps; ++n){

    Time.ResetStartTime();

    bc_disp = (double(n)+1.0)*delta_u;

    v.Import(lhs_u, *elasticityBVP->ImportToOverlapMap, Insert);
    w.Import(lhs_d, *phaseFieldBVP->ImportToOverlapMap, Insert);

    elasticityBVP->solve_u(matrix_u, rhs_u, lhs_u, v, w, *phaseFieldBVP->OverlapMap, ParametersList, bc_disp);

    elasticityBVP->updateDamageHistory(damageHistory, lhs_u);

    phaseFieldBVP->solve_d(matrix_d, rhs_d, lhs_d, ParametersList, damageHistory);

    if (Comm->MyPID()==0){
      std::cout << n << std::setw(15) << Time.ElapsedTime() << "\n";
    }

  }

  std::string dispfile = path + "displacement" + std::to_string(int(n_steps)) + ".mtx";
  std::string damgfile = path + "damage"       + std::to_string(int(n_steps)) + ".mtx";

  int error_u = elasticityBVP->print_solution(lhs_u, dispfile);
  int error_d = phaseFieldBVP->print_solution(lhs_d, damgfile);

}
