#ifndef DAMAGEPROBLEMREBAR_HPP
#define DAMAGEPROBLEMREBAR_HPP

#include "damageBVP.hpp"

class damageProblemRebar: public damageBVP
{

public:

  Epetra_SerialDenseVector gcs;

  damageProblemRebar(Epetra_Comm & comm, mesh & mesh_, double & gc, double & lc): damageBVP(comm, mesh_, lc){
    setup_dirichlet_conditions();
    gcs.Resize(2);
    gcs(0) = gc;
    gcs(1) = 100*gc;
  }
  ~damageProblemRebar(){
  }

  void get_fracture_energy(unsigned int & e_lid, unsigned int & gp, double & gc_){
    int e_gid = Mesh->local_cells[e_lid];
    int j     = Mesh->cells_physicalgroup[e_gid];
    gc_       = gcs(j);
  }

  void setup_dirichlet_conditions(){

  double tol = 1.0e-8;

  n_bc_dof = 0;
  double z;
  int node;
  for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
      node = Mesh->local_nodes[i];
      z    = Mesh->nodes_coord[3*node+2];
      if(z<=tol && z>=-tol){
          n_bc_dof++;
      }
      if(z<=1.0+tol && z>=1.0-tol){
          n_bc_dof++;
      }
  }

  int indbc = 0;
  dof_on_boundary = new int [n_bc_dof];
  for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
      node = Mesh->local_nodes[inode];
      z    = Mesh->nodes_coord[3*node+2];
      if(z<=tol && z>=-tol){
        dof_on_boundary[indbc] = inode;
        indbc++;
      }
      if(z<=1.0+tol && z>=1.0-tol){
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

  double tol = 1.0e-8;

  for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
      node = Mesh->local_nodes[inode];
      z    = Mesh->nodes_coord[3*node+2];
      if(z<=tol && z>=-tol){
          F[0][StandardMap->LID(node)] = 0.0;
      }
      if(z<=1.0+tol && z>=1.0-tol){
          F[0][StandardMap->LID(node)] = 0.0;
      }
  }
  ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
  }
};
#endif
