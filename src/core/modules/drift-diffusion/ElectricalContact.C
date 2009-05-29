// $Id$

#include "OhmicContact.h"
#include "SchottkyContact.h"
#include "LeakageCurrent.h"
#include "MaterialInterface.h"

#include "ElectricalContact.h"
#include "DriftDiffusionProperties.h"
#include "Boundary.h"
#include "Material.h"
#include "SimulationEnvironment.h"
#include "FowlerNordheim.h"


ElectricalContact*
ElectricalContact::create(const std::string& name,
    const ModelOptions& options)
{
  ElectricalContact* ct = NULL;

  if (name == "ohmic")
    ct = OhmicContact::create();
  else if (name == "schottky")
    ct = SchottkyContact::create();
  else if (name == "leakage_current")
    ct = LeakageCurrent::create();
  else if (name == "interface")
    ct = MaterialInterface::create();

  if (ct != NULL)
    ct->set_options(options);

  return ct;
}


void
ElectricalContact::do_init(void)
{
  if (get_option("zero_field", false))
    set_zero_derivative_bc(DriftDiffusionDefs::POTENTIAL);
  if (get_option("zero_grad_fermi_e", false))
    set_zero_derivative_bc(DriftDiffusionDefs::FERMIE);
  if (get_option("zero_grad_fermi_h", false))
    set_zero_derivative_bc(DriftDiffusionDefs::FERMIH);

  get_parameter("contact_resistance", _surfres);

  get_parameter("voltage", _voltage);

  _has_field_emission =
    get_option("calculate_field_emission", _has_field_emission);

  if (_has_field_emission)
  {
    _field_emission =
      new FowlerNordheim(get_option("work_function", 1.0));
    _real_contact = true;
  }

  determine_reference_material();

  // should be done in DriftDiffusionProperties
  get_reference_material().set_lattice_temperature(SimulationOptions::T);
  get_reference_material().calculate_equilibrium_properties();
  get_reference_material().setup_band_edges();
}


double
ElectricalContact::get_contact_voltage_drop(void) const
{
  double j = Constants::e * (get_normal_hole_flux() -
      get_normal_electron_flux());

  // a negative current means inflowing current
  return -_surfres * j;
}


void
ElectricalContact::determine_reference_material(void)
{
  Boundary* bd = get_boundary();
  const Device& dev = bd->get_environment()->get_device();
  const Material* mat;

  std::string ref = get_option("reference_material", "");
  if (!ref.empty())
  {
    mat = dev.get_material(ref);
    if (mat == NULL)
    {
      std::ostringstream os;
      os << "ElectricalContact: \'" << ref << "\' is not a valid"
        << " name for the reference material (region).";
      throw InitFailedException(os.str());
    }
    _reference_prop = dynamic_cast<DriftDiffusionProperties*>(
          mat->get_model(get_simulation_id()));
  }
  else
  {
    std::set<ID>::iterator it(bd->get_region_ids().begin());
    const std::set<ID>::iterator end(bd->get_region_ids().end());
    for ( ; it != end; ++it)
    {
      mat = dev.get_material(*it);
      DriftDiffusionProperties* prop =
        dynamic_cast<DriftDiffusionProperties*>(
            mat->get_model(get_simulation_id()));
      if (prop != NULL)
      {
        _reference_prop = prop;
        if (!prop->is_dielectric())
          break;
      }
    }
  }
}


double
ElectricalContact::calculate_field_emission(double F)
{
  double J = 0.0;
  if (_field_emission != NULL)
    J = _field_emission->get_emission_current(F);

  return J;
}
