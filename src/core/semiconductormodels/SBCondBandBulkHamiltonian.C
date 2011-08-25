// $Id$

#include "SBCondBandBulkHamiltonian.h"
#include "Semiconductor.h"
#include "Material.h"


SBCondBandBulkHamiltonian::~SBCondBandBulkHamiltonian(void)
{
}



void
SBCondBandBulkHamiltonian::prepare_submodels(void)
{
  assert(semiconductor == NULL);

  ModelOptions opt =  get_options();
  opt.delete_all_submodels();
  semiconductor = Semiconductor::create(get_material()->get_structure(), opt);
  add_submodel("semiconductor", semiconductor);
}



void
SBCondBandBulkHamiltonian::do_init(void)
{
  SBbulkHamiltonian::do_init();

  kp_bands.resize(1,0);

  kp_bands_map.insert(std::make_pair(0,0));

  calculate_for_init();

  get_option("model", "");
}




