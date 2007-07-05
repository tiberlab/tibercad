#include "MacrostrainBoundaryProperties.h"

 


MacrostrainBoundaryProperties::MacrostrainBoundaryProperties() : BoundaryProperties()
{


}

//==================================================================================//
MacrostrainBoundaryProperties *
MacrostrainBoundaryProperties::create(const std::string & name,  const ModelOptions &   options)
{
  MacrostrainBoundaryProperties* result = NULL;

  if (name == "substrate")
    result = MacrostrainSubstrate::create();
  else if (name == "pressure")
    result = MacrostrainPressure::create();

  if (result != NULL)
    result->set_options(options);

  return result;
}


//==================================================================================//

MacrostrainSubstrate::MacrostrainSubstrate()  :MacrostrainBoundaryProperties()
{
  material = NULL;
}


//===================================================================================//
MacrostrainSubstrate::~MacrostrainSubstrate()
{
  delete material;
}

//===================================================================================//
void MacrostrainSubstrate::do_init()
{
  const ModelOptions& options =	get_options ();

 

  const std::string name = options.get_option("material", "");

  if (name == "") throw InitFailedException("MacrostrainSubstrate: substrate material is not defined");


  material = Material::create (name, options);

  material->init();

  type = "substrate";

}


//====================================================================================//
MacrostrainPressure::MacrostrainPressure() :  MacrostrainBoundaryProperties()
{
  value = 0;

}

//====================================================================================//
void MacrostrainPressure::do_init()
{
  const ModelOptions& options =	get_options ();

  value = options.get_option("pressure", 0);

  type = "pressure";

}

