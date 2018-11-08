#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_RCP.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#include "ApplicationTemplate.hpp"

int main(int argc, char *argv[]){

  std::string xmlInFileName = "";

  Teuchos::CommandLineProcessor clp(false);
  clp.setOption("parameters",&xmlInFileName,"The XML file to read into a parameter list");
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

  Teuchos::RCP<Teuchos::ParameterList> Parameters = Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::updateParametersFromXmlFile(xmlInFileName, inoutArg(*Parameters));

  std::string filename = Teuchos::getParameter<std::string>(Parameters->sublist("Mesh"), "filename");
  TPF::mesh mesh(Comm, filename, 1.0);

  Teuchos::RCP<ApplicationTemplate> example = Teuchos::rcp(new ApplicationTemplate(Comm, mesh, *Parameters));

  example->setup_dirichlet_conditions();

  example->staggeredAlgorithmDirichletBC(*Parameters, true);

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
return 0;

}
