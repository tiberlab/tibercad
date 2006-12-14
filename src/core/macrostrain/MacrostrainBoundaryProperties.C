#include "MacrostrainBoundaryProperties.h"

MacrostrainBoundaryProperties::MacrostrainBoundaryProperties() : BoundaryProperties()
{


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

  const std::string name = options.get_option("material", "GaAs");

  material = Material::create (name, options);

  material->init();

  

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

}

