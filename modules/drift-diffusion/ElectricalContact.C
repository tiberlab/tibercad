// $Id$

#include "ElectricalContact.h"
#include "DriftDiffusionProperties.h"
#include "tibercad/physics/misc/ParticleDensity.h"
#include "tibercad/base/Initializer.h"


ElectricalContact::ElectricalContact(const ModelOptions& options)
  : DDInterfaceModel(options),
    _voltage(0.0),
    _surfres(0.0),
    _contact_fermilevel(0.0),
    _vrec_n(-1),
    _vrec_p(-1),
    _fixed_vrec_n(false),
    _fixed_vrec_p(false),
    _tunneling(false),
    _alpha_n(1.0),
    _phi_n(0.0),
    _mass_n(1.0),
    _alpha_p(1.0),
    _phi_p(0.0),
    _mass_p(1.0),
    _bias(0.0),
    _Vbi(0.0),
    _length(1.0)
{
  // this is a real contact
  has_current(true);
}


void
ElectricalContact::do_init(void)
{

  if (get_option("zero_field", false))
    set_type(0, NEUMANN);

  DDInterfaceModel::do_init();

  get_parameter("rec_velocity_n", _vrec_n);
  get_parameter("rec_velocity_p", _vrec_p);

  if (get_option("zero_grad_fermi_e", false) || (_vrec_n >= 0))
  {
    set_type(1, NEUMANN);
    _fixed_vrec_n = true;
  }
  else
    _vrec_n = 0.0;

  if (get_option("zero_grad_fermi_h", false) || (_vrec_p >= 0))
  {
    set_type(2, NEUMANN);
    _fixed_vrec_p = true;
  }
  else
    _vrec_p = 0.0;

  
  _tunneling = get_option("tunneling", false);
  if (_tunneling) 
  {
    get_parameter("alpha_n", _alpha_n);
    get_parameter("alpha_p", _alpha_p);
    get_parameter("mass_n", _mass_n);
    get_parameter("mass_p", _mass_p);
    get_parameter("phi_n", _phi_n);
    get_parameter("phi_p", _phi_p);
    get_parameter("bias", _bias);
    get_parameter("Vbi", _Vbi);
    get_parameter("length", _length);
  }
  

  get_parameter("contact_resistance", _surfres);

  get_parameter("voltage", _voltage);
}


void
ElectricalContact::do_compute(void)
{
  if (get_type(0) != NEUMANN)
    coeff_g(0) = _contact_fermilevel + get_inner_voltage();

  if (get_type(1) != NEUMANN)
    coeff_g(1) = get_inner_voltage();
  else
  {
    const PointData& pd = get_point_data();

    // we have to take the correct equilibrium density!
    BandProperties& cb = get_bulk_dd_properties()->get_conduction_band();
    cb.calculate(pd.electron_vt);
    std::pair<double, double> eldens(cb.get_density_and_derivative(0, _contact_fermilevel));
    double n0 = eldens.first;
    //ParticleDensity& el = get_dd_properties()->get_electrons();
    //el.set_classical_parameters(cb.get_effective_DOS(),
    //    get_dd_properties()->get_conduction_band_edge() - _contact_fermilevel, 0,
    //    pd.electron_vt);
    //double n0 = el.get_particle_density();

    //std::cout<<"vrec_n="<<_vrec_n<<" n="<<pd.electron_density<<" n0="<<n0<<std::endl;
    double Rn = _vrec_n * (pd.electron_density - n0);
    double dRn = _vrec_n * pd.electron_density_derivative;
    //std::cout<<"dRn = "<<dRn<<std::endl;
   
    if (_tunneling)
    {
      double Veff = _bias - _Vbi;
      if (Veff > 0.0 && _alpha_n > 0.0) 
      {
        //double A = Constants::e * Constants::e * 100 * 100 / (8 * M_PI * Constants::h);
        double B = 8 * M_PI * sqrt(2 * Constants::me * Constants::e) / (3 * Constants::h * 100);
        Rn -= _alpha_n * (Veff * Veff / (_length * _length)) * exp(-B * sqrt(_mass_n * _phi_n) * _phi_n  * _length/ Veff);
        //std::cout<<"B = "<<B<<" _alpha_n / _lenght^2 = "<<_alpha_n / (_length * _length)<<" b = "<<B * sqrt(_mass_n * _phi_n) * _phi_n  * _length<<" exp = "<<exp(-B * sqrt(_mass_n * _phi_n) * _phi_n  * _length/ Veff)<<" V^2 = "<<Veff * Veff<<std::endl;
      }
    }

    coeff_g(1) += Rn;
    jacobian(1, 0) += dRn;
    jacobian(1, 1) -= dRn;
  }

  if (get_type(2) != NEUMANN)
    coeff_g(2) = get_inner_voltage();
  else
  {
    const PointData& pd = get_point_data();

    // we have to take the correct equilibrium density!
    BandProperties& vb = get_bulk_dd_properties()->get_valence_band();
    vb.calculate(pd.hole_vt);
    std::pair<double, double> hldens(vb.get_density_and_derivative(0, _contact_fermilevel));
    double p0 = hldens.first;
    //ParticleDensity& hl = get_dd_properties()->get_holes();
    //vb.set_temperature(pd.hole_vt);
    //hl.set_classical_parameters(vb.get_effective_DOS(),
    //    _contact_fermilevel - get_dd_properties()->get_valence_band_edge(), 0,
    //    pd.hole_vt);
    //double p0 = hl.get_particle_density();

    //std::cout<<"vrec_p="<<_vrec_p<<" p="<<pd.hole_density<<" p0="<<p0<<std::endl;
    double Rp = _vrec_p * (pd.hole_density - p0);
    double dRp = _vrec_p * pd.hole_density_derivative;
    //std::cout<<"dRp ="<<dRp<<std::endl;

    if (_tunneling)
    {
      double Veff = _bias - _Vbi;
      if (Veff > 0.0 && _alpha_p > 0.0)
      {
        //double h_eV = Constants::h / Constants::e;
        //double A = Constants::e * Constants::e * 100 * 100 / (8 * M_PI * Constants::h);
        double B = 8 * M_PI * sqrt(2 * Constants::me * Constants::e) / (3 * Constants::h * 100);
        Rp -= _alpha_p * (Veff * Veff / (_length * _length)) * exp(-B * sqrt(_mass_p * _phi_p) * _phi_p * _length/ Veff);
      }
    }
 


    coeff_g(2) += Rp;
    jacobian(2, 0) += dRp;
    jacobian(2, 2) -= dRp;
  }

}

/*
double
ElectricalContact::get_contact_voltage_drop(void) const
{
  double j = Constants::e * (get_normal_hole_flux() -
      get_normal_electron_flux());

  // a negative current means inflowing current
  return -_surfres * j;
}
*/



