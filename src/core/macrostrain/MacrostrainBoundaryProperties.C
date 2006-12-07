#include "MacrostrainBoundaryProperties.h"

MacrostrainBoundaryProperties::MacrostrainBoundaryProperties() : BoundaryProperties()
{


}

//==================================================================================//

MacrostrainSubstrate::MacrostrainSubstrate()  :MacrostrainBoundaryProperties()
{
  crystal = NULL;
}


//===================================================================================//
MacrostrainSubstrate::~MacrostrainSubstrate()
{
  delete crystal;
}

//===================================================================================//
void MacrostrainSubstrate::do_init()
{
  const ModelOptions& options =	get_options ();

  const std::string name = options.get_option("material", "GaAs");

  // crystal = RotatedCrystal::create ( name, options ); 

  crystal->init();

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

