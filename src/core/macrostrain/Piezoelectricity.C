// $Id$

#include "Piezoelectricity.h"
#include "Material.h"
 
 
Piezoelectricity::Piezoelectricity(const ModelOptions& options)
 : PhysicalModelInterface(options)
{


}


Piezoelectricity* Piezoelectricity::create(const Material* mat, const ModelOptions& options)
{

  std::string structure = mat->get_structure();
  return dynamic_cast<Piezoelectricity*>(PhysicalModelInterface::create("piezo_" + structure,
      mat, options));

}
//-----------------------------------------------------------//

