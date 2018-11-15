/*
Brian Staber (brian.staber@gmail.com)
*/

#include "meshpp.hpp"

mesh::mesh(){
}

mesh::mesh(std::string & fileName_mesh, double scaling){
    read_gmsh(fileName_mesh, scaling);
}

mesh::mesh(Epetra_Comm & comm, std::string & fileName_mesh, double scaling){
    Comm = &comm;
    int MyPID = Comm->MyPID();
    int NumProc = Comm->NumProc();

    read_gmsh(fileName_mesh, scaling);

    epart = new idx_t[n_cells];
    npart = new idx_t[n_nodes];

    Comm->Barrier();
    if (Comm->NumProc()>1){
        if (MyPID==0){
            metis_part_mesh(NumProc);
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

    Comm->Broadcast(epart,n_cells,0);
    Comm->Broadcast(npart,n_nodes,0);
    Comm->Barrier();
    get_local_nodes(MyPID);
    get_cells_and_ghosts(MyPID);

    store_feinterp_cells();
}

mesh::~mesh(){
    delete[] epart;
    delete[] npart;
}

int mesh::read_gmsh(std::string & fileName_mesh, double scaling){
    int error = 0;
    std::ifstream meshfile;
    meshfile.open(fileName_mesh.c_str());

    if (!meshfile){
        std::cout << "*ERR* Can't open the input file\n ";
        error = 1;
        return error;
    }

    char buf[100],c;
    double xyz[3];
    unsigned int n_total_cells;
    unsigned int el_info;
    unsigned int num;
    unsigned int nbtag;
    unsigned int tag1;
    unsigned int tag2;
    unsigned int n_cells4 = 0;
    unsigned int node;

    unsigned int nodes_cells4[4];
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
        meshfile >> xyz[0];
        meshfile >> xyz[1];
        meshfile >> xyz[2];
        nodes_coord[3*i+0] = xyz[0]/scaling;
        nodes_coord[3*i+1] = xyz[1]/scaling;
        nodes_coord[3*i+2] = xyz[2]/scaling;
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
                cells_physicalgroup.push_back(tag1-1);
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

    n_cells4 = tetra4_nodes.size()/4;
    if (n_cells4==0){ // && n_cells10==0 && n_hexa8==0){
        std::cerr << "Your mesh is empty!\n";
    }
    else{
        n_cells = n_cells4;
        el_type = 4;
        cells_nodes.reserve(tetra4_nodes.size());
        cells_nodes = tetra4_nodes;

        gauss_points_tetra4(gauss_weight_cells,xi_cells,eta_cells,zeta_cells);
        n_gauss_cells = gauss_weight_cells.Length();
    }
    return error;
}

int mesh::read_boundary_file(std::string & fileName_bc, unsigned int & number_physical_groups){

    std::ifstream file_bc;
    file_bc.open(fileName_bc.c_str());
    if (!file_bc){
        std::cout << "*ERR* Can't open the input file called " << fileName_bc.c_str()  << "\n";
        return 1;
    }

    Epetra_IntSerialDenseVector input(number_physical_groups*n_nodes);
    for (unsigned int i=0; i<number_physical_groups*n_nodes; ++i){
        file_bc >> input(i);
    }

    nodes_to_boundaries.Reshape(n_local_nodes_without_ghosts,number_physical_groups);
    int node;
    for (unsigned int inode=0; inode<n_local_nodes_without_ghosts; ++inode){
        node = local_nodes[inode];
        for (unsigned pg=0; pg<number_physical_groups; ++pg){
        nodes_to_boundaries(inode,pg) = input(number_physical_groups*node+pg);
        }
    }
    file_bc.close();
    return 0;

}

int mesh::metis_part_mesh(int & NumProc){
    int check_PartMeshDual;

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

    idx_t *vwgt=NULL;
    idx_t *vsize=NULL;
    idx_t common = 3;
    idx_t nparts = NumProc;
    real_t *tpwgts=NULL;
    idx_t *options=NULL;
    idx_t objval;
    check_PartMeshDual = METIS_PartMeshDual(&ne, &nn, eptr, eind, vwgt, vsize, &common, &nparts, tpwgts, options, &objval, epart, npart);

    if (check_PartMeshDual==0){
        std::cout << "*ERR* An error occured with METIS.\n";
    }
    delete [] eptr;
    delete [] eind;
    delete [] vwgt;
    delete [] vsize;
    delete [] tpwgts;
    delete [] options;

    return check_PartMeshDual;
}

void mesh::get_local_nodes(int & MyPID){
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

void mesh::get_cells_and_ghosts(int & MyPID){

    Epetra_IntSerialDenseVector mynpart(Copy,npart,n_nodes);
    local_nodes = local_nodes_without_ghosts;
    local_dof = local_dof_without_ghosts;
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

Epetra_SerialDenseVector mesh::get_cartesian_coordinate(unsigned int & e_gid, unsigned int & gp){

  int node;
  int n_local_cells = n_local_cells;
  int n_gauss_cells = n_gauss_cells;
  double xi, eta, zeta;

  Epetra_SerialDenseVector vector_x(3);
  Epetra_SerialDenseVector shape_functions(el_type);
  Epetra_SerialDenseMatrix matrix_X(3,el_type);

  for (unsigned inode=0; inode<el_type; ++inode){
      node = cells_nodes[el_type*e_gid+inode];
      matrix_X(0,inode) = nodes_coord[3*node+0];
      matrix_X(1,inode) = nodes_coord[3*node+1];
      matrix_X(2,inode) = nodes_coord[3*node+2];
  }
  xi = xi_cells(gp); eta = eta_cells(gp); zeta = zeta_cells(gp);
  tetra4::shape_functions(shape_functions, xi, eta, zeta);
  vector_x.Multiply('N','N',1.0,matrix_X,shape_functions,0.0);
  return vector_x;
}

void mesh::store_feinterp_cells(){

    int node, eglob;
    double alpha, beta;
    Epetra_SerialDenseVector N(el_type);
    Epetra_SerialDenseMatrix JacobianMatrix(3,3), InverseJacobianMatrix(3,3), X(3,el_type), D(el_type,3), DX(el_type,3);

    local_rows.Resize(3*el_type*n_local_cells);
    vol_cells.Resize(n_local_cells);
    N_cells.Reshape(el_type,n_gauss_cells);
    detJac_cells.Reshape(n_local_cells,n_gauss_cells);
    DX_N_cells.Reshape(n_gauss_cells*el_type,n_local_cells);
    DY_N_cells.Reshape(n_gauss_cells*el_type,n_local_cells);
    DZ_N_cells.Reshape(n_gauss_cells*el_type,n_local_cells);

    for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
      tetra4::shape_functions(N, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
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
            local_rows(3*el_type*eloc+3*inode) = 3*node;
            local_rows(3*el_type*eloc+3*inode+1) = 3*node+1;
            local_rows(3*el_type*eloc+3*inode+2) = 3*node+2;
        }

        vol_cells(eloc) = 0.0;
        for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
            tetra4::d_shape_functions(D, xi_cells[gp], eta_cells[gp], zeta_cells[gp]);
            jacobian_matrix(X,D,JacobianMatrix);
            jacobian_det(JacobianMatrix,detJac_cells(eloc,gp));
            dX_shape_functions(D,JacobianMatrix,detJac_cells(eloc,gp),DX);
            vol_cells(eloc) += gauss_weight_cells(gp)*detJac_cells(eloc,gp);
            for (int inode=0; inode<el_type; ++inode){
                DX_N_cells(gp+n_gauss_cells*inode,eloc) = DX(inode,0);
                DY_N_cells(gp+n_gauss_cells*inode,eloc) = DX(inode,1);
                DZ_N_cells(gp+n_gauss_cells*inode,eloc) = DX(inode,2);
            }
        }

    }
}
