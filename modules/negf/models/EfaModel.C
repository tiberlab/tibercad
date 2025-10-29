// $Id$

#include "EfaModel.h"
#include "Messages.h"

#include "TiberModule.h"

EfaModel::EfaModel(const ModelOptions& options)
   : HamiltonianModel(options)
{
   //_inv_mass_crys(0);
   //_inv_mass(0);
}


void
EfaModel::do_init(void)
{
  _model = "efa";

  _simulation = get_option("simulation","none");

  if (_simulation=="none")
  {
    Messages::error("in hamiltonian submodel efa must define simulation");
    exit(1);
  }

  _num_bands = 0;

  _degeneracy = 1; 

}

