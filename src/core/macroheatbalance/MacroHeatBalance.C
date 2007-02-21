#include "MacroHeatBalance.h"
#include "DriftDiffusion.h"
#include "BoundaryProperties.h"
#include "PhysicalModel.h"

using namespace std;
//-----------------------------------------------------------------//


void MacroHeatBalance::parse_options( )
{ 

}

void MacroHeatBalance::do_init( ) 
{

}

void  MacroHeatBalance::do_solve()
{

}

MacroHeatBalance::~MacroHeatBalance()
{

}

MacroHeatBalance::MacroHeatBalance()
{
  _dd_simul = NULL;
}

PhysicalModel*   MacroHeatBalance :: create_physical_model (const ModelOptions &options) const 
                    throw (ModelErrorException)
{

}

BoundaryProperties* MacroHeatBalance :: create_boundary_model (const ModelOptions &options) const 
                    throw (ModelErrorException)

{

}

MacroHeatBalance*  MacroHeatBalance :: create (void)
{

}

void MacroHeatBalance :: build_nodal_results (const std::set< std::string > &variables, 
				     std::vector< double > &results, 
				     std::vector< std::string > &legend)
{

}

 
