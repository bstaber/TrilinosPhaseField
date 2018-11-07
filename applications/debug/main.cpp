#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "mesh.hpp"

int main(int argc, char *argv[]){

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif

  std::string filename = "/Users/brian/Documents/Github/TrilinosPhaseField/meshes/beam.msh";
  TPF::mesh Mesh(Comm, filename, 1.0);

  if (Comm.MyPID()==0){
    std::cout << "n_cells = " << Mesh.n_cells << "\n";
    std::cout << "n_nodes = " << Mesh.n_nodes << "\n";
    std::cout << "MyPID" << std::setw(15) << "n_local_cells" << std::setw(15) << "n_local_nodes" << std::setw(15) << "\n";
  }

  Comm.Barrier();
  std::cout << Comm.MyPID() << std::setw(15) << Mesh.n_local_cells << std::setw(15) << Mesh.n_local_nodes_without_ghosts << "\n";

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;

}
