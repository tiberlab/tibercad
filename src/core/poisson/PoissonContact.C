#include "PoissonContact.h"
#include "Dirichlet.h" 
#include "Neumann.h" 
//==================================================================================//
PoissonContact*
PoissonContact::create(const std::string & name,  const ModelOptions &   options)
{
  PoissonContact* result = NULL;



  if (name == "Dirichlet")
    result = Dirichlet::create(options);

  if (name == "Neumann")
    result = Neumann::create(options);
  
  return result;
}




