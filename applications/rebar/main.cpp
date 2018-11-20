/*
Brian Staber (brian.staber@gmail.com)
*/

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_RCP.hpp"
#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

#include "Teuchos_CommandLineProcessor.hpp"

#include "elasticProblemRebar.hpp"
#include "damageProblemRebar.hpp"
#include "staggeredAlgorithm.hpp"

int main(int argc, char *argv[]){

  std::string    xmlInFileName = "";

  Teuchos::CommandLineProcessor  clp(false);
  clp.setOption("xml-in-file",&xmlInFileName,"The XML file to read into a parameter list");
  clp.setDocString("TO DO.");

  Teuchos::CommandLineProcessor::EParseCommandLineReturn
  parse_return = clp.parse(argc,argv);
  if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) {
      std::cout << "\nEnd Result: TEST FAILED" << std::endl;
      return parse_return;
  }

#ifdef HAVE_MPI
MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Teuchos::RCP<Teuchos::ParameterList> paramList = Teuchos::rcp(new Teuchos::ParameterList);
  if(xmlInFileName.length()) {
      Teuchos::updateParametersFromXmlFile(xmlInFileName, inoutArg(*paramList));
  }
  if (Comm.MyPID()==0){
      paramList->print(std::cout,2,true,true);
  }

  std::string mesh_file = Teuchos::getParameter<std::string>(paramList->sublist("Mesh"), "mesh_file");
  double gc = Teuchos::getParameter<double>(paramList->sublist("Damage"), "gc");
  double lc = Teuchos::getParameter<double>(paramList->sublist("Damage"), "lc");

  mesh Mesh(Comm, mesh_file, 1.0);

  Teuchos::RCP<elasticProblemRebar> elasInterface   = Teuchos::rcp(new elasticProblemRebar(Comm, Mesh, *paramList));
  Teuchos::RCP<damageProblemRebar>  damageInterface = Teuchos::rcp(new damageProblemRebar(Comm, Mesh, gc, lc));
  Teuchos::RCP<staggeredAlgorithm> solver           = Teuchos::rcp(new staggeredAlgorithm(Comm, Mesh, *damageInterface, *elasInterface));

  solver->staggeredAlgorithmDirichletBC(*paramList, true);

  #ifdef HAVE_MPI
      MPI_Finalize();
  #endif
  return 0;
}
