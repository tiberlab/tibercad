// $Id$

#include "SimpleBand.h"
#include "Database.h"
#include "SimulationOptions.h"
#include "Messages.h"

#include "TiberModule.h"


SimpleBand::SimpleBand(const ModelOptions& options) :
  BandProperties(options),
  _reference_energy(0.0),
  _bandgap(0.0)
{
  _varshni[0] = 0.0;
  _varshni[1] = 0.0;
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
SimpleBand::read_database_alloy(void)
{
  const Database& db = get_database();

  db.set_section("bandgap");
  _bow_Eg = db.get("bow_Eg_G", 0.0);
}


void
SimpleBand::read_database(void)
{
  // when reading from the database, we use the same data
  // as for kp
  const Database& db = get_database();
  db.set_section("valenceband");
  _reference_energy = db.get("E_v", 0.0);

  if (get_option("particle", "el") == "el")
  {
    db.set_section("bandgap");
    _bandgap = db.get("Eg_G", 1e3);
    _varshni[0] = db.get("varshni_alpha_G", 0.0);
    _varshni[1] = db.get("varshni_beta_G", 0.0);

    double Eg_L = db.get("Eg_L", 1e3);
    double Eg_X = db.get("Eg_X", 1e3);

    db.set_section("conductionband");
    effective_mass() = db.get("m_G", 1.0);
    degeneracy() = 2;

    if (Eg_L < _bandgap)
    {
      _bandgap = Eg_L;

      db.set_section("bandgap");
      _varshni[0] = db.get("varshni_alpha_L", 0.0);
      _varshni[1] = db.get("varshni_beta_L", 0.0);

      db.set_section("conductionband");
      double m_L_t = db.get("m_L_t", 1.0);
      double m_L_l = db.get("m_L_l", 1.0);
      // 64 = 8 * 8, with 8 = valley degeneracy
      effective_mass() = std::pow(64 * m_L_t * m_L_t * m_L_l, 1.0/3.0 );
    }
    if (Eg_X < _bandgap)
    {
      _bandgap = Eg_X;

      db.set_section("bandgap");
      _varshni[0] = db.get("varshni_alpha_X", 0.0);
      _varshni[1] = db.get("varshni_beta_X", 0.0);

      db.set_section("conductionband");
      double m_X_t = db.get("m_X_t", 1.0);
      double m_X_l = db.get("m_X_l", 1.0);
      // 36 = 6 * 6, with 6 = valley degeneracy
      effective_mass() = std::pow(36 * m_X_t * m_X_t * m_X_l, 1.0/3.0 );
    }


  }
  else if (get_option("particle", "el") == "hl")
  {
    db.set_section("valenceband");
    degeneracy() = db.get("degeneracy", 4);
    effective_mass() = db.get("m_dos", 1.0);
  }
}

/*
void
SimpleBand::do_init_alloy(const PhysicalModelInterface* comp_A,
    const PhysicalModelInterface* comp_B, double xa)
{

  _modelA = dynamic_cast<const Semiconductor* >(comp_A);
  _modelB = dynamic_cast<const Semiconductor* >(comp_B);

  _xa = xa;
}
*/

void
SimpleBand::do_init(void)
{
  if (has_option("band_edge"))
  {
    _bandgap = 0.0;
    get_parameter("band_edge", _reference_energy);
    if (has_option("reference_energy"))
      Messages::warning("You defined both \'band_edge\' and \'reference_energy\'"
          " for semiconductor band parameters");

    if (has_option("band_gap"))
      Messages::warning("You defined both \'band_edge\' and \'band_gap\'"
          " for semiconductor band parameters");
  }
  get_parameter("reference_energy", _reference_energy);
  degeneracy() = 2;
  get_parameter("degeneracy", degeneracy());
  // add spin degeneracy to DOS mass
  double deg = std::pow(2.0, 2.0 / 3.0);

  if (get_option("particle", "el") == "el")
    get_parameter("band_gap", _bandgap);

  get_parameter("DOS_mass", effective_mass(), true,
      initializer(&SimpleBand::_set_mdos));
  _set_mdos(effective_mass());

  if (has_option("effective_DOS"))
  {
    get_parameter("effective_DOS", _eff_DOS, true,
        initializer(&SimpleBand::_set_mdos_from_Nc));
    _set_mdos_from_Nc(_eff_DOS);
  }

  do_calculate();
}




void
SimpleBand::do_calculate(void)
{
  band_edge() = _reference_energy;
  if (_bandgap > 0)
  {
    double T = get_lattice_temperature() / Constants::k_B;
    double gap = _bandgap;

    gap -= _varshni[0] * T * T / (T + _varshni[1]);

    band_edge() += gap;
  }
}
