/*
Brian Staber (brian.staber@gmail.com)
*/

#include "damageBVP.hpp"
#include "fepp.hpp"

damageBVP::damageBVP(Epetra_Comm & comm, mesh & mesh_, double & gc_, double & lc_):
gc(gc_), lc(lc_){

  Mesh               = &mesh_;
  Comm               = Mesh->Comm;

  StandardMap        = new Epetra_Map(-1, Mesh->n_local_nodes_without_ghosts, &Mesh->local_nodes_without_ghosts[0], 0, *Comm);
  OverlapMap         = new Epetra_Map(-1, Mesh->n_local_nodes, &Mesh->local_nodes[0], 0, *Comm);
  ImportToOverlapMap = new Epetra_Import(*OverlapMap, *StandardMap);

  create_FECrsGraph();
}

damageBVP::~damageBVP(){
}

void damageBVP::assemble(Epetra_FECrsMatrix & matrix, Epetra_FEVector & rhs,
                         Epetra_Vector & damageHistory, Epetra_Map & GaussMap){

  matrix.PutScalar(0.0);
  rhs.PutScalar(0.0);

  Epetra_SerialDenseVector shape_functions(Mesh->el_type);
  Epetra_SerialDenseVector fe(Mesh->el_type);

  Epetra_SerialDenseMatrix dx_shape_functions(3,Mesh->el_type);
  Epetra_SerialDenseMatrix ke(Mesh->el_type, Mesh->el_type);
  Epetra_SerialDenseMatrix me(Mesh->el_type, Mesh->el_type);

  double gauss_weight, hn, an;
  double bn = gc*lc;
  int eglob, id;
  int n_gauss_points = Mesh->n_gauss_cells;
  int * index = new int [Mesh->el_type];

  for (unsigned int eloc=0; eloc<Mesh->n_local_cells; ++eloc){
    eglob = Mesh->local_cells[eloc];
    for (int inode=0; inode<Mesh->el_type; ++inode){
      index[inode] = Mesh->cells_nodes[Mesh->el_type*eglob+inode];
      fe(inode) = 0.0;
      for (int jnode=0; jnode<Mesh->el_type; ++jnode){
        ke(inode,jnode) = 0.0;
        me(inode,jnode) = 0.0;
      }
    }

    for (unsigned int gp=0; gp<n_gauss_points; ++gp){
      gauss_weight = Mesh->gauss_weight_cells(gp);
      id = n_gauss_points*eglob+gp;
      hn = damageHistory[GaussMap.LID(id)];
      an = 2.0*hn + gc/double(lc);
      for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
          dx_shape_functions(0,inode) = Mesh->DX_N_cells(gp+n_gauss_points*inode,eloc);
          dx_shape_functions(1,inode) = Mesh->DY_N_cells(gp+n_gauss_points*inode,eloc);
          dx_shape_functions(2,inode) = Mesh->DZ_N_cells(gp+n_gauss_points*inode,eloc);
          shape_functions(inode) = Mesh->N_cells(inode,gp);
          fe(inode) += gauss_weight*2.0*hn*shape_functions(inode)*Mesh->detJac_cells(eloc,gp);
      }
      me.Multiply('N','T',an*gauss_weight*Mesh->detJac_cells(eloc,gp),shape_functions,shape_functions,1.0);
      ke.Multiply('T','N',bn*gauss_weight*Mesh->detJac_cells(eloc,gp),dx_shape_functions,dx_shape_functions,1.0);
    }
    ke += me;

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

void damageBVP::solve(Teuchos::ParameterList & Parameters,
                      Epetra_FECrsMatrix & matrix, Epetra_Vector & lhs, Epetra_FEVector & rhs,
                      Epetra_Vector & damageHistory,
                      Epetra_Map & GaussMap){

  double bc = 0.0;
  std::cout << "Step1\n";
  assemble(matrix, rhs, damageHistory, GaussMap);
  std::cout << "Step2\n";
  apply_dirichlet_conditions(matrix, rhs, bc);
  std::cout << "Step3\n";
  double tol   = Teuchos::getParameter<double>(Parameters.sublist("Aztec"), "AZ_tol");
  int max_iter = Teuchos::getParameter<int>(Parameters.sublist("Aztec"), "AZ_max_iter");

  lhs.PutScalar(0.0);

  Epetra_LinearProblem problem(&matrix, &lhs, &rhs);

  AztecOO solver(problem);
  solver.SetParameters(Parameters.sublist("Aztec"));
  solver.Iterate(max_iter, tol);
}

void damageBVP::create_FECrsGraph(){
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

int damageBVP::print_solution(Epetra_Vector & lhs, std::string fileName){
    int NumTargetElements = 0;
    if (Comm->MyPID()==0){
        NumTargetElements = Mesh->n_nodes;
    }
    Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
    Epetra_Export ExportOnRoot(*StandardMap,MapOnRoot);
    Epetra_MultiVector lhs_root(MapOnRoot,true);
    lhs_root.Export(lhs,ExportOnRoot,Insert);
    int error = EpetraExt::MultiVectorToMatrixMarketFile(fileName.c_str(),lhs_root,0,0,false);
    return error;
}
