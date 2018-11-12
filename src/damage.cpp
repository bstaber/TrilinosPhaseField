#include "damage.hpp"

TPF::damage::damage(Epetra_Comm & comm, mesh & mesh_, double & gc_, double & lc_):
gc(gc_), lc(lc_){

  Mesh = &mesh_;
  Comm = &comm;
}

TPF::damage::~damage(){
}

void TPF::damage::assemble(Epetra_FECrsMatrix & matrix, Epetra_FEVector & rhs,
                           Epetra_Vector & damageHistory){

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
    he = damageHistory[eloc];
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
        me.Multiply('N','T',gauss_weight,shape_functions,shape_functions,1.0);
      }
      me.Scale(ae*Mesh->detJac_cells(eloc));
      ke.Multiply('T','N',be*Mesh->vol_cells(eloc),dx_shape_functions,dx_shape_functions,0.0);
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

void TPF::damage::solve_d(Epetra_FECrsMatrix & A, Epetra_FEVector & rhs, Epetra_Vector & lhs,
                          Teuchos::ParameterList & Parameters, Epetra_Vector & damageHistory){

  assemble(A, rhs, damageHistory);

  double tol   = Teuchos::getParameter<double>(Parameters.sublist("Damage").sublist("Aztec"), "AZ_tol");
  int max_iter = Teuchos::getParameter<int>(Parameters.sublist("Damage").sublist("Aztec"), "AZ_max_iter");

  lhs.PutScalar(0.0);

  Epetra_LinearProblem problem(&A, &lhs, &rhs);
  AztecOO solver(problem);
  
  solver.SetParameters(Parameters.sublist("Damage").sublist("Aztec"));
  solver.Iterate(max_iter, tol);
}

int TPF::damage::print_solution(Epetra_Vector & lhs, std::string filename){
  int NumTargetElements = 0;
  if (Comm->MyPID()==0){
      NumTargetElements = Mesh->n_nodes;
  }
  Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
  Epetra_Export ExportOnRoot(*Mesh->StandardMapD,MapOnRoot);
  Epetra_MultiVector lhs_root(MapOnRoot,true);
  lhs_root.Export(lhs,ExportOnRoot,Insert);
  int error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
  return error;
}
