// $Id: DriftDiffusionProperties.C 4145 2015-10-02 11:53:20Z maufder $


#include "MasterEquationsProperties.h"
#include "SimulationInterface.h"
#include "Material.h"
#include "Database.h"

#include "Constants.h"
#include "Messages.h"
#include "TypeDefs.h"

#include "elem.h"


#include <cmath>


using namespace std;


MasterEquationsProperties::~MasterEquationsProperties(void)
{
}


MasterEquationsProperties::MasterEquationsProperties(const ModelOptions& options)
  : PhysicalModel(options),
    _elem(NULL),
    _statistics(TiberCad::BOLTZMANN),
    _conduction_band(NULL),
    _valence_band(NULL),
    _conduction_band_edge(NULL),
    _valence_band_edge(NULL)
    //_coupling(DriftDiffusionDefs::BOTH)
{

}



MasterEquationsProperties*
MasterEquationsProperties::create(const std::string& name, const Material* mat,
    const ModelOptions& options)
{
  //return PhysicalModel::create<MasterEquationsProperties>("me_bulk_" + name, mat, options);
}




void
MasterEquationsProperties::do_init(void)
{
  // hand the temperature interface over to the band parameter models
  _conduction_band->set_temperature_interface(_lattice_temp);
  _valence_band->set_temperature_interface(_lattice_temp);


  // set the temperature
  set_lattice_temperature(SimulationOptions::T);


}



void
MasterEquationsProperties::reinit(const Elem* elem)
{
  if (this->_elem != elem)
  {
    this->_elem = elem;

    lattice_vt = Constants::k_B * _lattice_temp.get_temperature(elem, elem->vertex_average());

    conduction_edge = _conduction_band_edge->get_band_edge();

    valence_edge = _valence_band_edge->get_band_edge();

    this->prepare_element_data();
  }
}



void
MasterEquationsProperties::set_conduction_band(MEBandProperties* cb)
{
  _conduction_band = cb;
}


void
MasterEquationsProperties::set_valence_band(MEBandProperties* vb)
{
  _valence_band = vb;
}


//void
//MasterEquationsProperties::get_conduction_band_edge(MEBandProperties* cb)
//{
//  _conduction_band_edge = cb;

//}


//void
//MasterEquationsProperties::get_valence_band_edge(MEBandProperties* vb)
//{
//  _valence_band_edge = vb;
//}
