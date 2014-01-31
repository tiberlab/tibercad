// $Id$

#include "EtbModel.h"
#include "Messages.h"

#include "TiberModule.h"

EtbModel::EtbModel(const ModelOptions& options)
   : HamiltonianModel(options)
{
   //_inv_mass_crys(0);
   //_inv_mass(0);
}


void
EtbModel::do_init(void)
{
  _model = "etb";

  _simulation = get_option("simulation","none");

  if (_simulation=="none")
  {
    Messages::error("in hamiltonian submodel etb must define simulation");
    exit(1);
  }

  _degeneracy = 1;
}

