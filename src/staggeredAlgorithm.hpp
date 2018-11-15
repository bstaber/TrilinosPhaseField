#ifndef STAGGEREDALGORITHM_HPP
#define STAGGEREDALGORITHM_HPP

#include "elasticBVP.hpp"
#include "damageBVP.hpp"

class staggeredAlgorithm
{

public:

  staggeredAlgorithm(damageBVP & damageInterface_, elasticBVP & elasticInterface_);
  ~staggeredAlgorithm();

  Epetra_Map constructGaussMap();

  void staggeredAlgorithmDirichletBC(Teuchos::ParameterList & ParametersList, bool print);

  damageBVP  * damageInterface;
  elasticBVP * elasticInterface;

  mesh        * Mesh;
  Epetra_Comm * Comm;

};
#endif
