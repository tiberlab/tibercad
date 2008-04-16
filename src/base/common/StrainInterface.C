// $Id: StrainInterface.C 827 2008-03-20 18:30:00Z gromano $

#include "StrainInterface.h"
#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "InitFailedException.h"

#include "elem.h"

#include "tensor.h"



StrainInterface::StrainInterface(void)
  : _simulation(NULL)
{
}



bool
StrainInterface::set_simulation(const std::string& name)
{
  bool answer = false;
  
  if (name != "")
  {
    _simulation = SimulationInterface::find_simulation(name);
    if (_simulation == NULL)
      throw InitFailedException("No such simulation found: " + name);

    // for now there are 9 variables from strain
    _strain_ids.resize(9);
    _id_set.clear();
    
    _strain_ids[0] = _simulation->get_variable_id("eps_xx");
    _strain_ids[1] = _simulation->get_variable_id("eps_yy");
    _strain_ids[2] = _simulation->get_variable_id("eps_zz");
    _strain_ids[3] = _simulation->get_variable_id("eps_xy");
    _strain_ids[4] = _simulation->get_variable_id("eps_yz");
    _strain_ids[5] = _simulation->get_variable_id("eps_xz");
    _strain_ids[6] = _simulation->get_variable_id("Px");
    _strain_ids[7] = _simulation->get_variable_id("Py");
    _strain_ids[8] = _simulation->get_variable_id("Pz");

    for (int i = 0; i < 9; i++)
      if (_strain_ids[i] == INVALID_ID)
        throw InitFailedException("Simulation " + name +
            " is missing strain related variables");
      else
        _id_set.insert(_strain_ids[i]);

    answer = true;
  }

  return answer;
}


void
StrainInterface::get_strain_data(const Elem* elem, Tensor2Sym& strain,
    RealVectorValue& polarization)
{
  assert(elem != NULL);

  strain = 0;
  polarization = 0;

  if (_simulation != NULL)
  {
    std::map<ID, double> data;
    bool ok = _simulation->get_solution(elem, elem->centroid(),
        _id_set, data);

    if (ok)
    {
      strain(1,1) = data[_strain_ids[0]];
      strain(2,2) = data[_strain_ids[1]];
      strain(3,3) = data[_strain_ids[2]];
      strain(2,1) = data[_strain_ids[3]];
      strain(3,2) = data[_strain_ids[4]];
      strain(3,1) = data[_strain_ids[5]];

      polarization(0) = data[_strain_ids[6]]; 
      polarization(1) = data[_strain_ids[7]]; 
      polarization(2) = data[_strain_ids[8]];
    }
  }
}


