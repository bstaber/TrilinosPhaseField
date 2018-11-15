/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef ELASTICBVP_HPP
#define ELASTICBVP_HPP

#include "linearizedElasticity.hpp"
#include "damageBVP.hpp"

class elasticBVP : public linearizedElasticity
{
public:
  Epetra_SerialDenseMatrix elasticity;
  Teuchos::RCP<damageBVP> damageInterface;

  double gc, lc;
  double E, nu, lambda, mu;

  elasticBVP(Epetra_Comm & comm, Teuchos::ParameterList & Parameters);
  ~elasticBVP();

  Epetra_Map constructGaussMap();

  void initialize(Epetra_Comm & comm, Teuchos::ParameterList & Parameters);

  void computeDisplacement(Teuchos::ParameterList & ParameterList,
                           Epetra_FECrsMatrix & matrix, Epetra_Vector & lhs, Epetra_FEVector & rhs,
                           Epetra_Vector & phi, Epetra_Map & OverlapMapD, double & bc_disp);

  void updateDamageHistory(Epetra_Vector & damageHistory,
                           Epetra_Vector & displacement,
                           Epetra_Map & GaussMap);

  void staggeredAlgorithmDirichletBC(Teuchos::ParameterList & ParametersList,
                                     bool print);

  void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, double & phi_e, Epetra_SerialDenseVector & epsilon, Epetra_SerialDenseMatrix & tangent_matrix);

};
#endif
