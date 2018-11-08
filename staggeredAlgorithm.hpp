#ifndef STAGGEREDALGORITHM_HPP
#define STAGGEREDALGORITHM_HPP

#include "elasticity.hpp"
#include "damage.hpp"

class staggeredAlgorithm : public elasticity
{

public:
  Epetra_SerialDenseMatrix elasticity_matrix;
  Teuchos::RCP<damage> damageInterface;

  double gc, lc;
  double E, nu, lambda, mu;

  Epetra_Vector * damageSolutionOverlaped;

  staggeredAlgorithm();
  ~staggeredAlgorithm();

  void initialize(Epetra_Comm & comm, Teuchos::ParameterList & Parameters);

  void computeDisplacement(Teuchos::ParameterList & ParameterList,
                           Epetra_FECrsMatrix & matrix, Epetra_Vector & lhs, Epetra_FEVector & rhs,
                           double & bc_disp);

  void updateDamageHistory(Epetra_Vector & damageHistory,
                           Epetra_Vector & displacement);

  void staggeredAlgorithmDirichletBC(Teuchos::ParameterList & ParametersList, bool print);

  void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, Epetra_SerialDenseMatrix & tangent_matrix);

};
#endif
