// $Id$

#include "SimpleBand.h"
#include "Database.h"
#include "SimulationOptions.h"


TIBER_MODULE(SimpleBand, band_properties, simple)


SimpleBand::SimpleBand(const ModelOptions& options) :
  BandProperties(options)
{

}



SimpleBand::~SimpleBand(void)
{

}


void
SimpleBand::_set_mdos(double& mdos)
{
  double deg = std::pow(2.0, 2.0 / 3.0);
  effective_mass() = deg * mdos;
}



void
SimpleBand::_set_mdos_from_Nc(double& Nc)
{
  double kT = Constants::k_B * SimulationOptions::T;
  effective_mass() =
      std::pow(Nc / get_dos_factor(), 2.0/3.0) / kT;
}


void
SimpleBand::read_database(void)
{
  Database& db = get_database();
  db.set_section("band_properties/simple");

  if (get_option("particle", "el") == "el")
  {
    band_edge() = db.get("Ec", 0.0, true);
    degeneracy() = db.get("g_el", 2, true);
    effective_mass() = db.get("DOS_mass_el", -1.0);
    if (effective_mass() < 0)
    {
      _eff_DOS = db.get("effective_DOS_el", 1e20, true);
      _set_mdos_from_Nc(_eff_DOS);
    }
    else
    _set_mdos(effective_mass());

  }
  else if (get_option("particle", "el") == "hl")
  {
    band_edge() = db.get("Ev", 0.0, true);
    degeneracy() = db.get("g_hl", 2, true);
    effective_mass() = db.get("DOS_mass_hl", -1.0);
    if (effective_mass() < 0)
    {
      _eff_DOS = db.get("effective_DOS_hl", 1e20, true);
      _set_mdos_from_Nc(_eff_DOS);
    }
    else
    _set_mdos(effective_mass());

  }
}



void
SimpleBand::do_init(void)
{
  get_parameter("band_edge", band_edge());
  get_parameter("reference_energy", band_edge());
  degeneracy() = 2;
  get_parameter("degeneracy", degeneracy());
  // add spin degeneracy to DOS mass
  double deg = std::pow(2.0, 2.0 / 3.0);
  effective_mass() = 1;
  get_parameter("DOS_mass", effective_mass(), true,
      initializer(&SimpleBand::_set_mdos));
  _set_mdos(effective_mass());

  if (has_option("effective_DOS"))
  {
    _eff_DOS = 1e20;
    get_parameter("effective_DOS", _eff_DOS, true,
        initializer(&SimpleBand::_set_mdos_from_Nc));
    _set_mdos_from_Nc(_eff_DOS);
  }
}
