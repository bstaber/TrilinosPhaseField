#ifndef APPLICATIONTEMPLATE_HPP
#define APPLICATIONTEMPLATE_HPP

#include "staggeredAlgorithm.hpp"

class ApplicationTemplate : public TPF::staggeredAlgorithm
{

public:
  ApplicationTemplate(Epetra_Comm & comm, TPF::mesh & mesh_, Teuchos::ParameterList & Parameters):
  staggeredAlgorithm(comm, mesh_){

  }
  ~ApplicationTemplate(){

  }

  void setup_dirichlet_conditions(){

  }

  void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){

  }

  void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, Epetra_SerialDenseMatrix & tangent_matrix){

  }

};
#endif
