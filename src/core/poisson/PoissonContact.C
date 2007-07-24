#include "PoissonContact.h"
#include "Dirichlet.h" 
#include "Neumann.h" 
//==================================================================================//
PoissonContact*
PoissonContact::create(const std::string & name,  const ModelOptions &   options)
{
  PoissonContact* result = NULL;



  if (name == "Dirichlet")
    result = Dirichlet::create();

  if (name == "Neumann")
    result = Neumann::create();
  


  if (result != NULL)
    result->set_options(options);

  return result;
}




