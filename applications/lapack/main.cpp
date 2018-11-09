#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_LAPACK.h"
#include "Epetra_SerialDenseVector.h"

int main(int argc, char *argv[]){

  #ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
  #else
    Epetra_SerialComm Comm;
  #endif

  Epetra_LAPACK Lapack;
  Epetra_SerialDenseVector A(6);
  A(0) = 0.8147; A(1) = 0.9096; A(3) = 0.2027;
                 A(2) = 0.6324; A(4) = 0.3222;
                                A(5) = 0.9575;

  char jobz = 'V';
  char uplo = 'U';

  double * ap   = new double[6]; // packed lower or upper triangular part of the matrix
  double * w    = new double[3]; // eigenvalues
  double * z    = new double[9]; // eigenvectors

  Epetra_SerialDenseVector eigv(3);
  Epetra_SerialDenseVector eigw(9);

  double * WORK = new double[6];
  int    INFO;

  ap = &A[0];

  Lapack.SPEV(jobz, uplo, 3, ap, &eigv[0], &eigw[0], 3, WORK, &INFO);

  std::cout << "Eigv = " << eigv << "\n";
  std::cout << "Eigw = " << eigw << "\n";

  #ifdef HAVE_MPI
      MPI_Finalize();
  #endif
  return 0;

}
