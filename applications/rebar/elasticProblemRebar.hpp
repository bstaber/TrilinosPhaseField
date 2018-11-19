/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef ELASTICPROBLEMREBAR_HPP
#define ELASTICPROBLEMREBAR_HPP

#include "elasticBVP.hpp"

class elasticProblemRebar : public elasticBVP
{
public:

  Epetra_SerialDenseMatrix elasticity;
  Epetra_SerialDenseVector lambda, mu, E, nu;

  elasticProblemRebar(Epetra_Comm & comm, mesh & mesh_, Teuchos::ParameterList & Parameters): elasticBVP(comm, mesh_, Parameters){
    setup_dirichlet_conditions();

    E.Resize(2);
    nu.Resize(2);
    lambda.Resize(2);
    mu.Resize(2);

    E(0)  = Teuchos::getParameter<double>(Parameters.sublist("Elasticity"), "young");
    nu(0) = Teuchos::getParameter<double>(Parameters.sublist("Elasticity"), "poisson");

    E(1)  = 10.0*E(0);
    nu(1) = 10.0*nu(0);

    for (unsigned int j=0; j<2; ++j){
      lambda(j) = E(j)*nu(j)/((1.0+nu(j))*(1.0-2.0*nu(j)));
      mu(j)     = E(j)/(2.0*(1.0+nu(j)));
    }

    elasticity.Reshape(6,6);
  }

  ~elasticProblemRebar(){
  }

  void get_elasticity_tensor(unsigned int & e_lid, unsigned int & gp, double & phi_e, Epetra_SerialDenseVector & epsilon, Epetra_SerialDenseMatrix & tangent_matrix){

    int e_gid = Mesh->local_cells[e_lid];
    int j     = Mesh->cells_physicalgroup[e_gid];

    double c11 = lambda(j)+2.0*mu(j);
    double c12 = lambda(j);
    double c44 = mu(j);

    elasticity(0,0) = c11; elasticity(0,1) = c12; elasticity(0,2) = c12; elasticity(0,3) = 0.0; elasticity(0,4) = 0.0; elasticity(0,5) = 0.0;
    elasticity(1,0) = c12; elasticity(1,1) = c11; elasticity(1,2) = c12; elasticity(1,3) = 0.0; elasticity(1,4) = 0.0; elasticity(1,5) = 0.0;
    elasticity(2,0) = c12; elasticity(2,1) = c12; elasticity(2,2) = c11; elasticity(2,3) = 0.0; elasticity(2,4) = 0.0; elasticity(2,5) = 0.0;
    elasticity(3,0) = 0.0; elasticity(3,1) = 0.0; elasticity(3,2) = 0.0; elasticity(3,3) = c44; elasticity(3,4) = 0.0; elasticity(3,5) = 0.0;
    elasticity(4,0) = 0.0; elasticity(4,1) = 0.0; elasticity(4,2) = 0.0; elasticity(4,3) = 0.0; elasticity(4,4) = c44; elasticity(4,5) = 0.0;
    elasticity(5,0) = 0.0; elasticity(5,1) = 0.0; elasticity(5,2) = 0.0; elasticity(5,3) = 0.0; elasticity(5,4) = 0.0; elasticity(5,5) = c44;

    double g = (1.0-phi_e)*(1.0-phi_e) + 1.0e-6;
    tangent_matrix = elasticity;
    tangent_matrix.Scale(g);

  }

  void updateDamageHistory(Epetra_Vector & damageHistory, Epetra_Vector & displacement, Epetra_Map & GaussMap){

    Epetra_Vector u(*OverlapMap);
    u.Import(displacement, *ImportToOverlapMap, Insert);

    int n_gauss_points = Mesh->n_gauss_cells;

    double trepsilon, trepsilon2, potential, history;

    Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
    Epetra_SerialDenseMatrix matrix_B(6,3*Mesh->el_type);
    Epetra_SerialDenseVector cells_u(3*Mesh->el_type);
    Epetra_SerialDenseVector epsilon(6);

    int egid, node, id, phase;
    for (unsigned int elid=0; elid<Mesh->n_local_cells; ++elid){
      egid = Mesh->local_cells[elid];
      phase = Mesh->cells_physicalgroup[egid];
      for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
        node = Mesh->cells_nodes[Mesh->el_type*egid+inode];
        for (unsigned int ddl=0; ddl<3; ++ddl){
          id = 3*node+ddl;
          cells_u(3*inode+ddl) = u[OverlapMap->LID(id)];
        }
      }
      for (unsigned int gp=0; gp<n_gauss_points; ++gp){
          for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
              dx_shape_functions(inode,0) = Mesh->DX_N_cells(gp+n_gauss_points*inode,elid);
              dx_shape_functions(inode,1) = Mesh->DY_N_cells(gp+n_gauss_points*inode,elid);
              dx_shape_functions(inode,2) = Mesh->DZ_N_cells(gp+n_gauss_points*inode,elid);
          }

          compute_B_matrices(dx_shape_functions,matrix_B);
          epsilon.Multiply('N','N',1.0,matrix_B,cells_u,0.0);

          trepsilon  = epsilon(0) + epsilon(1) + epsilon(2);
          trepsilon2 = epsilon(0)*epsilon(0) + epsilon(1)*epsilon(1) + epsilon(2)*epsilon(2) +
                       0.5*epsilon(3)*epsilon(3) + 0.5*epsilon(4)*epsilon(4) + 0.5*epsilon(5)*epsilon(5);

          id        = n_gauss_points*egid+gp;
          history   = damageHistory[GaussMap.LID(id)];
          potential = (lambda(phase)/2.0)*trepsilon*trepsilon + mu(phase)*trepsilon2;
          if (potential>history){
            damageHistory[GaussMap.LID(id)] = potential;
          }
      }
    }

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
            n_bc_dof+=3;
        }
        if(z<=1.0+tol && z>=1.0-tol){
            n_bc_dof+=1;
        }
    }

    int indbc = 0;
    dof_on_boundary = new int [n_bc_dof];
    for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
        node = Mesh->local_nodes[inode];
        z    = Mesh->nodes_coord[3*node+2];
        if (z<=tol && z>=-tol){
            dof_on_boundary[indbc+0] = 3*inode+0;
            dof_on_boundary[indbc+1] = 3*inode+1;
            dof_on_boundary[indbc+2] = 3*inode+2;
            indbc+=3;
        }
        if (z<=1.0+tol && z>=1.0-tol){
            dof_on_boundary[indbc+0] = 3*inode+2;
            indbc+=1;
        }
    }
  }

  void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double & displacement){

    double tol = 1.0e-8;

    Epetra_MultiVector v(*StandardMap,true);
    v.PutScalar(0.0);

    int node;
    double z;
    for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
        node = Mesh->local_nodes[inode];
        z    = Mesh->nodes_coord[3*node+2];
        if (z<=1.0+tol && z>=1.0-tol){
            v[0][StandardMap->LID(3*node+2)] = displacement;
        }
    }
    Epetra_MultiVector rhs_dir(*StandardMap,true);
    K.Apply(v,rhs_dir);
    F.Update(-1.0,rhs_dir,1.0);

    for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
        node = Mesh->local_nodes[inode];
        z    = Mesh->nodes_coord[3*node+2];
        if (z<=tol && z>=-tol){
            F[0][StandardMap->LID(3*node+0)] = 0.0;
            F[0][StandardMap->LID(3*node+1)] = 0.0;
            F[0][StandardMap->LID(3*node+2)] = 0.0;
        }
        if (z<=1.0+tol && z>=1.0-tol){
            F[0][StandardMap->LID(3*node+2)] = displacement;
        }
    }
    ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
  }
};
#endif
