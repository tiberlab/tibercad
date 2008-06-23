// $Id$


#include "ExcitonProperties.h"
#include "Constants.h"

#include "elem.h"


ExcitonProperties::~ExcitonProperties(void)
{
}


ExcitonProperties::ExcitonProperties(void)
  : _elem(NULL),
    _statistics(TiberCad::BOLTZMANN)
{
}


ExcitonProperties*
ExcitonProperties::create(const std::string& name,
    const ModelOptions& options)
{
  return dynamic_cast<ExcitonProperties*>(
      PhysicalModelInterface::create("exmodel_" + name, options));
}



void
ExcitonProperties::reinit(const Elem* elem)
{
  if (this->_elem != elem)
  {
    this->_elem = elem;
    lattice_vt = Constants::k_B *
      _lattice_temp.get_temperature(elem, elem->centroid());
    
    this->prepare_element_data();
  }
}

