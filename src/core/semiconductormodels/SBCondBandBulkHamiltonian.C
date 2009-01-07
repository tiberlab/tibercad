// $Id$

#include "SBCondBandBulkHamiltonian.h"
#include "Semiconductor.h"
#include "Material.h"


SBCondBandBulkHamiltonian::~SBCondBandBulkHamiltonian()
{
  PhysicalModelInterface::destroy(semiconductor);
}






void
SBCondBandBulkHamiltonian::do_init(void)
{
  SBbulkHamiltonian::do_init();

  kp_bands.resize(1,0);

  kp_bands_map.insert(std::make_pair(0,0));

  const ModelOptions& opt =  get_options ();

  PhysicalModelInterface::destroy(semiconductor);
  semiconductor = Semiconductor::create(get_material()->get_structure(), opt);
  semiconductor->set_material(get_material());
  semiconductor->init();

}




void
SBCondBandBulkHamiltonian::do_init_alloy(const PhysicalModelInterface *comp_A,
    const PhysicalModelInterface *comp_B, double xa)
{
  const SBCondBandBulkHamiltonian* matA =
    dynamic_cast< const SBCondBandBulkHamiltonian*> (comp_A);

  const SBCondBandBulkHamiltonian* matB =
    dynamic_cast< const SBCondBandBulkHamiltonian*> (comp_B);

  PhysicalModelInterface::destroy(semiconductor);
  semiconductor = static_cast<Semiconductor*>(matA->semiconductor->copy());
  assert(semiconductor != NULL);
  semiconductor->set_material(get_material());

  //Added by G.Penazzi
 kp_bands.resize(1,0);

  kp_bands_map.insert(std::make_pair(0,0));
  //--------------------

  semiconductor->init_alloy(matA->semiconductor, matB->semiconductor, xa);

  calculate_for_init();

}

