#include "mesh.hpp"

TPF::mesh::mesh(Epetra_Comm & comm, std::string & filename, double scaling){

    Comm        = &comm;
    int MyPID   = Comm->MyPID();
    int NumProc = Comm->NumProc();

    read_gmsh(filename, scaling);

    partition(MyPID, NumProc);

    set_gauss_points(4, gauss_weight_cells, xi_cells, eta_cells, zeta_cells);

    computeFE();

    build_EpetraMaps();
    build_FECrsGraphs();
}

TPF::mesh::~mesh(){
    delete[] epart;
    delete[] npart;
}

int TPF::mesh::read_gmsh(std::string & filename, double scaling){

    std::ifstream meshfile;
    meshfile.open(filename.c_str());

    if (!meshfile){
        std::cout << "*ERR* Can't open the input file\n ";
        return 1;
    }

    char   buf[100], c;
    double coord;

    unsigned int n_total_cells;
    unsigned int el_info;
    unsigned int num;
    unsigned int nbtag;
    unsigned int tag1;
    unsigned int tag2;

    unsigned int n_cells4 = 0;
    unsigned int node;

    std::vector<int> tetra4_nodes;

    meshfile.getline(buf,100);
    meshfile.getline(buf,100);
    meshfile.getline(buf,100);
    meshfile.getline(buf,100);
    meshfile >> n_nodes;
    meshfile.get(c);

    nodes_coord.reserve(3*n_nodes);

    for (int i=0; i<n_nodes; ++i){
        meshfile >> num;
        meshfile >> coord;
        nodes_coord[3*i+0] = coord/scaling;
        meshfile >> coord;
        nodes_coord[3*i+1] = coord/scaling;
        meshfile >> coord;
        nodes_coord[3*i+2] = coord/scaling;
    }

    meshfile.getline(buf,100);
    meshfile.getline(buf,100);
    meshfile.getline(buf,100);
    meshfile >> n_total_cells;
    meshfile.getline(buf,100);

    for (int i=0; i<n_total_cells; ++i){

        meshfile >> num;
        meshfile >> el_info;
        meshfile >> nbtag;
        meshfile >> tag1;
        meshfile >> tag2;

        switch (el_info) {
            case 1:
                for (unsigned int inode=0; inode<2; ++inode){
                    meshfile >> node;
                }
                break;
            case 2:
                for (unsigned int inode=0; inode<3; ++inode){
                    meshfile >> node;
                }
                break;
            case 3:
                for (unsigned int inode=0; inode<4; ++inode){
                    meshfile >> node;
                }
                break;
            case 4:
                for (unsigned int inode=0; inode<4; ++inode){
                    meshfile >> node;
                    tetra4_nodes.push_back(node-1);
                }
                break;
            case 5:
                for (unsigned int inode=0; inode<8; ++inode){
                    meshfile >> node;
                }
                break;
            case 9:
                for (unsigned int inode=0; inode<6; ++inode){
                    meshfile >> node;
                }
                break;
            case 11:
                for (unsigned int inode=0; inode<10; ++inode){
                    meshfile >> node;
                }
                break;
            case 15:
                for (unsigned int inode=0; inode<1; ++inode){
                    meshfile >> node;
                }
                break;
            default:
                std::cout << "Element not supported encountered: el_info = " << el_info << "\n";
                break;
        };
    }
    meshfile.close();

    n_cells = tetra4_nodes.size()/4;
    if (n_cells>0){
      cells_nodes.reserve(tetra4_nodes.size());
      cells_nodes = tetra4_nodes;
    }
    else{
      std::cout << "Mesh is empty.\n";
    }
    return 0;
}

void TPF::mesh::partition(int & MyPID, int & NumProc){

  epart = new idx_t[n_cells];
  npart = new idx_t[n_nodes];

  Comm->Barrier();
  if (Comm->NumProc()>1){
      if (MyPID==0){
          int flag = metis(NumProc);
      }
  }
  else{
      for (unsigned int e=0; e<n_cells; ++e){
          epart[e] = 0;
      }
      for (unsigned int n=0; n<n_nodes; ++n){
          npart[n] = 0;
      }
  }

  Comm->Broadcast(epart, n_cells, 0);
  Comm->Broadcast(npart, n_nodes, 0);
  Comm->Barrier();

  get_local_nodes(MyPID);
  get_cells_and_ghosts(MyPID);
}

int TPF::mesh::metis(int & NumProc){

    idx_t ne = n_cells;
    idx_t nn = n_nodes;
    idx_t *eptr, *eind;
    eptr = new idx_t[n_cells+1];
    eind = new idx_t[el_type*n_cells];

    for (unsigned int i=0; i<n_cells; ++i){
        for (unsigned int inode=0; inode<el_type; ++inode){
            eind[el_type*i+inode] = cells_nodes[el_type*i+inode];
        }
        eptr[i] = el_type*i;
    }
    eptr[n_cells] = el_type*n_cells;

    idx_t *vwgt    = NULL;
    idx_t *vsize   = NULL;
    idx_t common   = 3;
    idx_t nparts   = NumProc;
    real_t *tpwgts = NULL;
    idx_t *options = NULL;
    idx_t objval;

    int flag = METIS_PartMeshDual(&ne, &nn, eptr, eind, vwgt, vsize, &common, &nparts, tpwgts, options, &objval, epart, npart);

    if (flag==0){
        std::cout << "*ERR* An error occured with METIS.\n";
    }

    delete [] eptr;
    delete [] eind;
    delete [] vwgt;
    delete [] vsize;
    delete [] tpwgts;
    delete [] options;

    return flag;
}

void TPF::mesh::get_local_nodes(int & MyPID){
    for (unsigned int j=0; j<n_nodes; ++j){
        if (npart[j]==MyPID){
            local_nodes_without_ghosts.push_back(j);
            local_dof_without_ghosts.push_back(3*j+0);
            local_dof_without_ghosts.push_back(3*j+1);
            local_dof_without_ghosts.push_back(3*j+2);
        }
    }
    n_local_nodes_without_ghosts = local_nodes_without_ghosts.size();
}

void TPF::mesh::get_cells_and_ghosts(int & MyPID){

    Epetra_IntSerialDenseVector mynpart(Copy,npart,n_nodes);
    local_nodes = local_nodes_without_ghosts;
    local_dof   = local_dof_without_ghosts;
    int node;

    for (unsigned int i=0; i<n_cells; ++i){
        if (epart[i]==MyPID){
            local_cells.push_back(i);
            for (unsigned inode=0; inode<el_type; ++inode){
                node = cells_nodes[el_type*i+inode];
                if (mynpart[node]!=MyPID){
                    mynpart[node]=MyPID;
                    local_nodes.push_back(node);
                    local_dof.push_back(3*node+0);
                    local_dof.push_back(3*node+1);
                    local_dof.push_back(3*node+2);
                }
            }
        }
    }

    n_local_nodes = local_nodes.size();
    n_local_cells = local_cells.size();
}


void TPF::mesh::shape_functions(Epetra_SerialDenseVector & N, double & xi, double & eta, double & zeta){
    double lambda = 1.0 - xi - eta - zeta;
    N(0) = lambda;
    N(1) = xi;
    N(2) = eta;
    N(3) = zeta;
}

void TPF::mesh::d_shape_functions(Epetra_SerialDenseMatrix & D){
    D(0,0) = -1.0; D(0,1) = -1.0; D(0,2) = -1.0;
    D(1,0) = 1.0;  D(1,1) = 0.0;  D(1,2) = 0.0;
    D(2,0) = 0.0;  D(2,1) = 1.0;  D(2,2) = 0.0;
    D(3,0) = 0.0;  D(3,1) = 0.0;  D(3,2) = 1.0;
}

void TPF::mesh::dX_shape_functions(Epetra_SerialDenseMatrix & D, Epetra_SerialDenseMatrix JacobianMatrix, double & jac, Epetra_SerialDenseMatrix & DX)
{
    Epetra_SerialDenseMatrix InverseJacobianMatrix(3,3);

    InverseJacobianMatrix(0,0) = (1.0/jac)*(JacobianMatrix(1,1)*JacobianMatrix(2,2)-JacobianMatrix(1,2)*JacobianMatrix(2,1));
    InverseJacobianMatrix(0,1) = (1.0/jac)*(JacobianMatrix(0,2)*JacobianMatrix(2,1)-JacobianMatrix(0,1)*JacobianMatrix(2,2));
    InverseJacobianMatrix(0,2) = (1.0/jac)*(JacobianMatrix(0,1)*JacobianMatrix(1,2)-JacobianMatrix(0,2)*JacobianMatrix(1,1));

    InverseJacobianMatrix(1,0) = (1.0/jac)*(JacobianMatrix(1,2)*JacobianMatrix(2,0)-JacobianMatrix(1,0)*JacobianMatrix(2,2));
    InverseJacobianMatrix(1,1) = (1.0/jac)*(JacobianMatrix(0,0)*JacobianMatrix(2,2)-JacobianMatrix(0,2)*JacobianMatrix(2,0));
    InverseJacobianMatrix(1,2) = (1.0/jac)*(JacobianMatrix(0,2)*JacobianMatrix(1,0)-JacobianMatrix(0,0)*JacobianMatrix(1,2));

    InverseJacobianMatrix(2,0) = (1.0/jac)*(JacobianMatrix(1,0)*JacobianMatrix(2,1)-JacobianMatrix(1,1)*JacobianMatrix(2,0));
    InverseJacobianMatrix(2,1) = (1.0/jac)*(JacobianMatrix(0,1)*JacobianMatrix(2,0)-JacobianMatrix(0,0)*JacobianMatrix(2,1));
    InverseJacobianMatrix(2,2) = (1.0/jac)*(JacobianMatrix(0,0)*JacobianMatrix(1,1)-JacobianMatrix(0,1)*JacobianMatrix(1,0));

    DX.Multiply('N','N',1.0,D,InverseJacobianMatrix,0.0);
}

void TPF::mesh::jacobian_matrix(Epetra_SerialDenseMatrix & X, Epetra_SerialDenseMatrix & D, Epetra_SerialDenseMatrix & JacobianMatrix){
  JacobianMatrix.Multiply('N','N',1.0,X,D,0.0);
}

void TPF::mesh::jacobian_det(Epetra_SerialDenseMatrix & JacobianMatrix, double & jac){
    jac = std::fabs(JacobianMatrix(0,0)*JacobianMatrix(1,1)*JacobianMatrix(2,2)-JacobianMatrix(0,0)*JacobianMatrix(1,2)*JacobianMatrix(2,1)-JacobianMatrix(0,1)*JacobianMatrix(1,0)*JacobianMatrix(2,2)+JacobianMatrix(0,1)*JacobianMatrix(1,2)*JacobianMatrix(2,0)+JacobianMatrix(0,2)*JacobianMatrix(1,0)*JacobianMatrix(2,1)-JacobianMatrix(0,2)*JacobianMatrix(1,1)*JacobianMatrix(2,0) );
}

void TPF::mesh::set_gauss_points(int ngp, Epetra_SerialDenseVector & weight, Epetra_SerialDenseVector & xi, Epetra_SerialDenseVector & eta, Epetra_SerialDenseVector & zeta){
    double alpha = (5.0 - sqrt(5.0))/20.0;
    double beta  = (5.0 + 3.0*sqrt(5.0))/20.0;
    switch (ngp){
      case 1:
        weight.Resize(1);
        xi.Resize(1); eta.Resize(1); zeta.Resize(1);
        weight(0) = 1.0/6.0;
        xi(0) = 1.0/4.0;
        eta(0) = 1.0/4.0;
        zeta(0) = 1.0/4.0;
        break;
      case 4:
        weight.Resize(4); xi.Resize(4); eta.Resize(4); zeta.Resize(4);
        weight(0) = 1.0/24.0; weight(1) = 1.0/24.0; weight(2) = 1.0/24.0; weight(3) = 1.0/24.0;
        xi(0) = alpha; eta(0) = alpha; zeta(0) = alpha;
        xi(1) = alpha; eta(1) = alpha; zeta(1) = beta;
        xi(2) = alpha; eta(2) = beta;  zeta(2) = alpha;
        xi(3) = beta;  eta(3) = alpha; zeta(3) = alpha;
        break;
      case 5:
        weight.Resize(5); xi.Resize(5); eta.Resize(5); zeta.Resize(5);
        xi(0) = 0.25; xi(1) = 0.50; xi(3) = 0.1666666666666667; xi(4) = 0.1666666666666667; xi(5) = 0.1666666666666667;
        eta(0) = 0.25; eta(1) = 0.1666666666666667; eta(2) = 0.1666666666666667; eta(3) = 0.1666666666666667; eta(4) = 0.50;
        zeta(0) = 0.25; zeta(1) = 0.1666666666666667; zeta(2) = 0.1666666666666667; zeta(3) = 0.50; zeta(4) = 0.1666666666666667;
        weight(0) = -0.8000000000000000/6.0; weight(1) = 0.45/6.0; weight(2) = 0.45/6.0; weight(3) = 0.45/6.0; weight(4) =0.45/6.0;
        break;
      case 11:
        weight.Resize(11); xi.Resize(11); eta.Resize(11); zeta.Resize(11);
        xi(0) = 0.25; xi(1) = 0.7857142857142857; xi(2) = 0.0714285714285714; xi(3) = 0.0714285714285714; xi(4) = 0.0714285714285714; xi(5) = 0.1005964238332008; xi(6) = 0.3994035761667992; xi(7) = 0.3994035761667992; xi(8) = 0.3994035761667992; xi(9) = 0.1005964238332008; xi(10) = 0.1005964238332008;
        eta(0) = 0.25; eta(1) = 0.0714285714285714; eta(2) = 0.0714285714285714; eta(3) = 0.0714285714285714; eta(4) = 0.7857142857142857; eta(5) = 0.3994035761667992; eta(6) = 0.1005964238332008; eta(7) = 0.3994035761667992; eta(8) = 0.1005964238332008; eta(9) = 0.3994035761667992; eta(10) = 0.1005964238332008;
        zeta(0) = 0.25; zeta(1) = 0.0714285714285714; zeta(2) = 0.0714285714285714; zeta(3) = 0.7857142857142857; zeta(4) = 0.0714285714285714; zeta(5) = 0.3994035761667992; zeta(6) = 0.3994035761667992; zeta(7) = 0.1005964238332008; zeta(8) = 0.1005964238332008; zeta(9) = 0.1005964238332008; zeta(10) = 0.3994035761667992;
        weight(0) = -0.0789333333333333/6.0; weight(1) = 0.0457333333333333/6.0; weight(2) = 0.0457333333333333/6.0; weight(3) = 0.0457333333333333/6.0; weight(4) = 0.0457333333333333/6.0; weight(5) = 0.1493333333333333/6.0; weight(6) = 0.1493333333333333/6.0; weight(7) = 0.1493333333333333/6.0; weight(8) = 0.1493333333333333/6.0; weight(9) = 0.1493333333333333/6.0; weight(10) = 0.1493333333333333/6.0;
        break;
      default:
        std::cout << "Please pick a number of gauss points equal to 1, 4, 5, or 11: 4 is the default.\n";
        weight.Resize(4); xi.Resize(4); eta.Resize(4); zeta.Resize(4);
        weight(0) = 1.0/24.0; weight(1) = 1.0/24.0; weight(2) = 1.0/24.0; weight(3) = 1.0/24.0;
        xi(0) = alpha; eta(0) = alpha; zeta(0) = alpha;
        xi(1) = alpha; eta(1) = alpha; zeta(1) = beta;
        xi(2) = alpha; eta(2) = beta;  zeta(2) = alpha;
        xi(3) = beta;  eta(3) = alpha; zeta(3) = alpha;
        break;
    };
}

void TPF::mesh::computeFE(){

    int node, eglob;
    double alpha, beta;
    Epetra_SerialDenseVector N(el_type);
    Epetra_SerialDenseMatrix JacobianMatrix(3,3), InverseJacobianMatrix(3,3), X(3,el_type), D(el_type,3), DX(el_type,3);

    local_rows.Resize(3*el_type*n_local_cells);
    vol_cells.Resize(n_local_cells);
    detJac_cells.Resize(n_local_cells);

    N_cells.Reshape(el_type,n_gauss_cells);
    DX_N_cells.Reshape(el_type,n_local_cells);
    DY_N_cells.Reshape(el_type,n_local_cells);
    DZ_N_cells.Reshape(el_type,n_local_cells);

    d_shape_functions(D);

    for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
      shape_functions(N, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
      for (int inode=0; inode<el_type; ++inode){
        N_cells(inode,gp) = N(inode);
      }
    }

    for (unsigned int eloc=0; eloc<n_local_cells; ++eloc){

        eglob = local_cells[eloc];

        for (unsigned int inode=0; inode<el_type; ++inode){
            node = cells_nodes[el_type*eglob+inode];
            X(0,inode) = nodes_coord[3*node];
            X(1,inode) = nodes_coord[3*node+1];
            X(2,inode) = nodes_coord[3*node+2];
            local_rows(3*el_type*eloc+3*inode)   = 3*node;
            local_rows(3*el_type*eloc+3*inode+1) = 3*node+1;
            local_rows(3*el_type*eloc+3*inode+2) = 3*node+2;
        }

        jacobian_matrix(X,D,JacobianMatrix);
        jacobian_det(JacobianMatrix,detJac_cells(eloc));
        dX_shape_functions(D,JacobianMatrix,detJac_cells(eloc),DX);
        vol_cells(eloc) = (1.0/6.0)*detJac_cells(eloc);

        for (int inode=0; inode<el_type; ++inode){
          DX_N_cells(inode,eloc) = DX(inode,0);
          DY_N_cells(inode,eloc) = DX(inode,1);
          DZ_N_cells(inode,eloc) = DX(inode,2);
        }
    }
}

void TPF::mesh::build_EpetraMaps(){

  StandardMapU        = new Epetra_Map(-1, 3*n_local_nodes_without_ghosts, &local_dof_without_ghosts[0], 0, *Comm);
  OverlapMapU         = new Epetra_Map(-1, 3*n_local_nodes, &local_dof[0], 0, *Comm);
  ImportToOverlapMapU = new Epetra_Import(*OverlapMapU, *StandardMapU);

  StandardMapD        = new Epetra_Map(-1, n_local_nodes_without_ghosts, &local_nodes_without_ghosts[0], 0, *Comm);
  OverlapMapD         = new Epetra_Map(-1, n_local_nodes, &local_nodes[0], 0, *Comm);
  ImportToOverlapMapD = new Epetra_Import(*OverlapMapD, *StandardMapD);

}

void TPF::mesh::build_FECrsGraphs(){

  FEGraphU = new Epetra_FECrsGraph(Copy,*StandardMapU,100);
  int * indexU = new int [3*el_type];

  int eglob, node;
  for (int e_lid=0; e_lid<n_local_cells; ++e_lid){
      eglob = local_cells[e_lid];
      for (int inode=0; inode<el_type; ++inode){
          node = cells_nodes[el_type*eglob+inode];
          for (int ddl=0; ddl<3; ++ddl){
              indexU[3*inode+ddl] = 3*node+ddl;
          }
      }
      for (int i=0; i<3*el_type; ++i){
          for (int j=0; j<3*el_type; ++j){
              FEGraphU->InsertGlobalIndices(1, &indexU[i], 1, &indexU[j]);
          }
      }

  }
  Comm->Barrier();
  FEGraphU->GlobalAssemble();
  delete [] indexU;

  FEGraphD = new Epetra_FECrsGraph(Copy,*StandardMapD,100);
  int * indexD = new int [el_type];

  for (int eloc=0; eloc<n_local_cells; ++eloc){
      eglob = local_cells[eloc];
      for (int inode=0; inode<el_type; ++inode){
          node = cells_nodes[el_type*eglob+inode];
          indexD[inode] = node;
      }

      for (int i=0; i<el_type; ++i){
          for (int j=0; j<el_type; ++j){
              FEGraphD->InsertGlobalIndices(1, &indexD[i], 1, &indexD[j]);
          }
      }
  }
  Comm->Barrier();
  FEGraphD->GlobalAssemble();
  delete [] indexD;

}
