#include "damage.hpp"

TPF::damage::damage(Epetra_Comm & comm, mesh & mesh_, double & gc_, double & lc_):
gc(gc_), lc(lc_){

  Mesh = &mesh_;
  Comm = &comm;

  StandardMap = new Epetra_Map(-1, Mesh->n_local_nodes_without_ghosts, &Mesh->local_nodes_without_ghosts[0], 0, *Comm);
  OverlapMap  = new Epetra_Map(-1, Mesh->n_local_nodes, &Mesh->local_nodes[0], 0, *Comm);
  ImportToOverlapMap = new Epetra_Import(*OverlapMap, *StandardMap);

  create_FECrsGraph();
}

TPF::damage::~damage(){
}

void TPF::damage::assemble(Epetra_FECrsMatrix & matrix, Epetra_FEVector & rhs,
                           Epetra_MultiVector & damageHistory){

  matrix.PutScalar(0.0);
  rhs.PutScalar(0.0);

  Epetra_SerialDenseVector shape_functions(Mesh->el_type);
  Epetra_SerialDenseVector fe(Mesh->el_type);

  Epetra_SerialDenseMatrix dx_shape_functions(3,Mesh->el_type);
  Epetra_SerialDenseMatrix ke(Mesh->el_type, Mesh->el_type);
  Epetra_SerialDenseMatrix me(Mesh->el_type, Mesh->el_type);

  int eglob, id;
  int n_gauss_points = Mesh->n_gauss_cells;
  int * index = new int [Mesh->el_type];

  double gauss_weight;
  double be = gc*lc;
  double ae, he;

  for (unsigned int eloc=0; eloc<Mesh->n_local_cells; ++eloc){

    eglob = Mesh->local_cells[eloc];
    he = damageHistory[0][eloc];
    ae = 2.0*he + gc/double(lc);

    for (int inode=0; inode<Mesh->el_type; ++inode){
      index[inode] = Mesh->cells_nodes[Mesh->el_type*eglob+inode];
      dx_shape_functions(0,inode) = Mesh->DX_N_cells(inode,eloc);
      dx_shape_functions(1,inode) = Mesh->DY_N_cells(inode,eloc);
      dx_shape_functions(2,inode) = Mesh->DZ_N_cells(inode,eloc);
      fe(inode) = 0.0;
      for (int jnode; jnode<Mesh->el_type; ++jnode){
        ke(inode,jnode) = 0.0;
        me(inode,jnode) = 0.0;
      }

      for (unsigned int gp=0; gp<n_gauss_points; ++gp){
        gauss_weight = Mesh->gauss_weight_cells(gp);
        for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
          shape_functions(inode) = Mesh->N_cells(inode,gp);
          fe(inode) += gauss_weight*2.0*he*shape_functions(inode)*Mesh->detJac_cells(eloc);
        }
        me.Multiply('N','T',ae*gauss_weight*Mesh->detJac_cells(eloc),shape_functions,shape_functions,1.0);
      }

      ke.Multiply('T','N',be*gauss_weight*Mesh->detJac_cells(eloc),dx_shape_functions,dx_shape_functions,0.0);
      ke += me;
    }

    for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
      rhs.SumIntoGlobalValues(1, &index[inode], &fe(inode));
      for (unsigned int jnode=0; jnode<Mesh->el_type; ++jnode){
        matrix.SumIntoGlobalValues(1, &index[inode], 1, &index[jnode], &ke(inode,jnode));
      }
    }

  }

  Comm->Barrier();

  matrix.GlobalAssemble();
  matrix.FillComplete();
  rhs.GlobalAssemble();

  delete [] index;

}

void TPF::damage::solve(Teuchos::ParameterList & Parameters,
                        Epetra_FECrsMatrix & matrix, Epetra_Vector & lhs, Epetra_FEVector & rhs,
                        Epetra_Vector & damageHistory){

  assemble(matrix, rhs, damageHistory);

  double tol   = Teuchos::getParameter<double>(Parameters.sublist("Aztec"), "AZ_tol");
  int max_iter = Teuchos::getParameter<int>(Parameters.sublist("Aztec"), "AZ_max_iter");

  lhs.PutScalar(0.0);

  Epetra_LinearProblem problem(&matrix, &lhs, &rhs);

  AztecOO solver(problem);
  solver.SetParameters(Parameters.sublist("Aztec"));
  solver.Iterate(max_iter, tol);
}

void TPF::damage::create_FECrsGraph(){

  FEGraph = new Epetra_FECrsGraph(Copy,*StandardMap,100);
  int eglob, node;
  int *index;
  index = new int [Mesh->el_type];

  for (int eloc=0; eloc<Mesh->n_local_cells; ++eloc){
      eglob = Mesh->local_cells[eloc];
      for (int inode=0; inode<Mesh->el_type; ++inode){
          node = Mesh->cells_nodes[Mesh->el_type*eglob+inode];
          index[inode] = node;
      }

      for (int i=0; i<Mesh->el_type; ++i){
          for (int j=0; j<Mesh->el_type; ++j){
              FEGraph->InsertGlobalIndices(1, &index[i], 1, &index[j]);
          }
      }
  }
  Comm->Barrier();
  FEGraph->GlobalAssemble();
  delete[] index;
}
