#ifndef DAMAGEPROBLEM_HPP
#define DAMAGEPROBLEM_HPP

#include "damageBVP.hpp"

class damageProblem: public damageBVP
{

public:

  damageProblem(Epetra_Comm & comm, mesh & mesh_, double & gc_, double & lc_): damageBVP(comm, mesh_, gc_, lc_){
  }
  ~damageProblem(){
  }

  void get_fracture_energy(unsigned int & e_lid, unsigned int & gp, double & gc){

  }

  void setup_dirichlet_conditions(){
    //std::cout << "No essential boundary conditions.\n";
  }

  void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){
    //std::cout << "No essential boundary conditions.\n";
  }

};
#endif
