#ifndef BLOCK_SINGLE_NOTCH_HPP
#define BLOCK_SINGLE_NOTCH_HPP

#include "staggeredAlgorithm.hpp"

class block_single_notch : public TPF::staggeredAlgorithm
{

public:

  unsigned int n_bc_dof;
  int * dof_on_boundary;

  block_single_notch(Epetra_Comm & comm, TPF::mesh & mesh_, Teuchos::ParameterList & Parameters):
  staggeredAlgorithm(comm, mesh_){
    // add steps if needed
  }
  ~block_single_notch(){
    // delete any pointers that you introduced in this file
  }

  void setup_dirichlet_conditions(){

    // example of how the setup should be implemented
    // this should be simplified or hidden in next versions

    n_bc_dof = 0;
    double z;
    int node;

    // get the total number of prescribed dof
    for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
        node = Mesh->local_nodes[i];
        z    = Mesh->nodes_coord[3*node+2]; // where do you want
        if(z<=1.0e-6 && z>=-1.0e-6){        // to apply Dirichlet conditions ?
            n_bc_dof+=3;                    // how many dof at this location ?
        }
        if(z<=10+1.0e-6 && z>=10-1.0e-6){   // how many dof at this location ?
            n_bc_dof+=1;
        }
    }

    // store the prescibred dof into a vector dof_on_boundary
    int indbc = 0;
    dof_on_boundary = new int [n_bc_dof];
    for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
        node = Mesh->local_nodes[inode];
        z    = Mesh->nodes_coord[3*node+2];
        if (z<=1.0e-6 && z>=-1.0e-6){
            dof_on_boundary[indbc+0] = 3*inode+0;
            dof_on_boundary[indbc+1] = 3*inode+1;
            dof_on_boundary[indbc+2] = 3*inode+2;
            indbc+=3;
        }
        if (z<=10+1.0e-6 && z>=10-1.0e-6){
            dof_on_boundary[indbc+0] = 3*inode+2;
            indbc+=1;
        }
    }
  }

  void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){

    // example of how dirichlet conditions should be applied
    // this should be simplified or hidden in next versions

    Epetra_MultiVector v(*Mesh->StandardMapU,true);
    v.PutScalar(0.0);

    int node;
    double z;
    for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
        node = Mesh->local_nodes[inode];
        z    = Mesh->nodes_coord[3*node+2];
        if (z<=10+1.0e-6 && z>=10-1.0e-6){
            v[0][Mesh->StandardMapU->LID(3*node+2)] = displacement;
        }
    }

    Epetra_MultiVector rhs_dir(*Mesh->StandardMapU,true);
    K.Apply(v,rhs_dir);
    F.Update(-1.0,rhs_dir,1.0);

    for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
        node = Mesh->local_nodes[inode];
        z    = Mesh->nodes_coord[3*node+2];
        if (z<=1.0e-6 && z>=-1.0e-6){
            F[0][Mesh->StandardMapU->LID(3*node+0)] = 0.0;
            F[0][Mesh->StandardMapU->LID(3*node+1)] = 0.0;
            F[0][Mesh->StandardMapU->LID(3*node+2)] = 0.0;
        }
        if (z<=10+1.0e-6 && z>=10-1.0e-6){
            F[0][Mesh->StandardMapU->LID(3*node+2)] = displacement;
        }
    }
    //}
    ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
  }

};
#endif
