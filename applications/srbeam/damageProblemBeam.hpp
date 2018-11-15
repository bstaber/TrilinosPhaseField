#ifndef DAMAGEPROBLEMBEAM_HPP
#define DAMAGEPROBLEMBEAM_HPP

#include "damageBVP.hpp"

class damageProblemBeam: public damageBVP
{

public:

  damageProblemBeam(Epetra_Comm & comm, mesh & mesh_, double & gc_, double & lc_): damageBVP(comm, mesh_, gc_, lc_){
    setup_dirichlet_conditions();
  }
  ~damageProblemBeam(){
  }

  void setup_dirichlet_conditions(){

  n_bc_dof = 0;
  double z;
  int node;
  for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
      node = Mesh->local_nodes[i];
      z    = Mesh->nodes_coord[3*node+2];
      if(z==0.0 || z==1.0){
          n_bc_dof++;
      }
  }

  int indbc = 0;
  dof_on_boundary = new int [n_bc_dof];
  for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
      node = Mesh->local_nodes[inode];
      z    = Mesh->nodes_coord[3*node+2];
      if (z==0.0 || z==1.0){
          dof_on_boundary[indbc] = inode;
          indbc++;
      }
  }

  }

  void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){

  Epetra_MultiVector v(*StandardMap,true);
  v.PutScalar(0.0);

  int node;
  double z;

  Epetra_MultiVector rhs_dir(*StandardMap,true);
  K.Apply(v,rhs_dir);
  F.Update(-1.0,rhs_dir,1.0);

  for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
      node = Mesh->local_nodes[inode];
      z    = Mesh->nodes_coord[3*node+2];
      if (z==0.0 || z==1.0){
          F[0][StandardMap->LID(node)] = 0.0;
      }
  }
  ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
  }
};
#endif
