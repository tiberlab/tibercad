// $Id$

#ifndef _SEMICONDUCTORMODEL_H_
#define _SEMICONDUCTORMODEL_H_

#include "SimulationOptions.h"
#include "Dopant.h"
#include "DriftDiffusionProperties.h"
#include "Constants.h"

#include <vector>
#include <string>
#include <cmath>

extern "C" {
  #include "petscsnes.h"
}

// forward declarations
//class Point;
class Elem;
class DDsemiconductor;

//! A generic semiconductor model
/*!
 * This model implements Boltzmann and Fermi-Dirac statistics,
 * constant mobility
 */
class SemiconductorModel : public DriftDiffusionProperties
{
  public:

    SemiconductorModel(void);
    SemiconductorModel(const SemiconductorModel& model);
    virtual ~SemiconductorModel(void);

    //! Set n-type doping
    void set_n_dopant(const Dopant& dopant);

    //! Set p-type doping
    void set_p_dopant(const Dopant& dopant);

    void set_mobilities(double mu_e, double mu_h)
    {
      _conduction_band.mobility = mu_e;
      _valence_band.mobility = mu_h;
    };
    
    //! \copydoc DriftDiffusionProperties::calculate_equilibrium_properties()
    virtual void calculate_equilibrium_properties(
        int coupling = DriftDiffusionDefs::BOTH,
        double temperature = SimulationOptions::T);

    //! Add recombination model to be used.
    void add_recombination_model(DriftDiffusionDefs::RecombinationModel
        recomb_model);

    //! Remove a recombination model.
    void remove_recombination_model(DriftDiffusionDefs::RecombinationModel
        recomb_model);

    //! Set Shockley-Read-Hall recombination parameters
    void set_SRH_parameters(double tau_n, double tau_p);

    //! Set exciton generation rate parameter
    void set_exciton_generation_rate_parameter(double g)
      { _exciton_gen_param = g; };

    //! Set Direct Recombination parameters
    void set_direct_rec_parameters(double C)
      { _direct_rec_param = C; };

    //! Set the coupling type
    void set_coupling_type(DriftDiffusionDefs::Coupling coupling)
      { _coupling = coupling; };

    void set_data_file(const std::string& filename)
      { _filename = filename; };

    virtual void read_database(const Dummy&);

    //! \deprecated { Create parameters for an alloy }
    /*!
     * \deprecated { This method will live as long as the database is
     * not used yet.}
     */
    virtual void build_alloy(const std::string& component2,
        const std::string& bowing_params, double content);

    /*! \copydoc DriftDiffusionProperties::calculate_all()
     * 
     * This implementation models a simple semiconductor
     */
    virtual void calculate_all(double potential,
      double fermi_e, double fermi_h, const Point& p);
    
    /*!
     * @returns the donor density
     */
    double get_donor_density(void) const;

    /*!
     * @returns the acceptor density
     */
    double get_acceptor_density(void) const;

    /*!
     * @returns the band gap
     */
    double get_band_gap(void) const;

    //! Get the conduction band edge
    /*!
     * \return the conduction band edge
     */
    double get_conduction_band_edge(void) const;

    //! Get the valence band edge
    /*!
     * \return the valence band edge
     */
    double get_valence_band_edge(void) const;

    //! Get the exciton generation rate
    double get_exciton_generation_rate(void) const;


    void print_info(void) const;
    
  protected:
    
    /*!
     * This structure holds the basic properties of a band for given
     * conditions (temp etc.)
     */
    struct BandProperties
    {
      //! The effective mass for the DOS
      /*!
       * It includes any degeneration, i.e. also spin
       */
      double effective_mass;
      //! The effective density of states
      double effective_DOS;
      //! The band edge energy
      double band_edge;
      //! The low field mobility
      double mobility;
    };

    //! Get the physical semiconductor model
    /*!
     * Derived classes will need to access the physical model, e.g. to
     * set the strain.
     */
    DDsemiconductor* get_physical_model(void)
      { return _bulk_model; };

    //! Get the conduction band properties
    BandProperties& get_conduction_band_properties(void)
      { return _conduction_band; };

    //! Get the valence band properties
    BandProperties& get_valence_band_properties(void)
      { return _valence_band; };

    //! Extract the band properties from _bulk_model
    /*!
     * This method looks for the band extrema and puts the effective
     * mass, band edges etc. into the BandProperties structure
     */
    void extract_band_properties(void);

    //! \copydoc DriftDiffusionProperties::prepare_element_data()
    virtual void prepare_element_data(void);

    //! Setup the band edge data
    /*!
     * This implementation calculates the effective density of states
     * and sets the band edges.
     */
    void setup_band_edges(void);

    //! Get the constant factor to calculate the effective density of states
    /*!
     * \return the factor 2 * pow(2 * PI / h^2)^1.5
     */
    double get_DOS_factor(void) const;

    //! Calculate ionized donor and acceptor densities and derivatives
    /*!
     * \param[in] arg_e = (Ef_e - Ec + potential) / kT
     * \param[in] arg_h = (Ev - Ef_h - potential) / kT
     * \param[in] kT the thermal voltage
     * \param[out] Nd the ionized donor density
     * \param[out] dNd the derivative of \p Nd
     */
    void calculate_ionized_donors(double arg_e, double kT,
        double& Nd, double& dNd);

    void calculate_ionized_acceptors(double arg_h, double kT,
        double& Na, double& dNa);

    //! Calculate Shockley-Read-Hall recombination
    void calculate_SRH_recombination(void);

    //! Calculate Auger recombination
    void calculate_Auger_recombination(void);

    //! Calculate direct recombination
    void calculate_direct_recombination(void);

    //! Calculate exciton generation
    void calculate_exciton_generation(void);

    //! Calculates all properties according to the type of coupling
    /*!
     * In this simple model we assume no dependence on the coordinates
     */
    template <int coupling>
    void calculate_all(double potential, double fermi_e, double fermi_h);

    //! Set the object to unprepared state
    void set_to_unprepared(void);

  private:

    //! A flag to tell the state of this object
    /*!
     * \c true means that all data is prepared and ready for use
     */
    bool _is_prepared;

    //! The physical model for this semiconductor
    /*!
     * The physical model is based on an effective mass approximation
     */
    DDsemiconductor* _bulk_model;

    /*!
     * The factor 2 * pow(2 * PI / h^2)^1.5
     * for calculating the effective density of states
     */
    double _DOS_factor;


    std::string _filename;

    
    //! The recombination models used
    int _recombination;


    //! Type of coupling (particles) we want to study
    /*!
     * This can be one of \c ELECTRONS, \c HOLES or \c BOTH
     */
    int _coupling;


    /**
     * The band properties
     */
    BandProperties _conduction_band;
    BandProperties _valence_band;

    Dopant _n_dopant;
    Dopant _p_dopant;

    double _electron_recombination_time;
    double _hole_recombination_time;

    double _direct_rec_param;

    double _exciton_gen_param;


    static PetscErrorCode jacobian(SNES snes, Vec x,
        Mat *jac, Mat *B, MatStructure *flag, void *sc);

    static PetscErrorCode function(SNES snes, Vec x, Vec f, void *sc);
};


//
// inline member functions
//

inline
void
SemiconductorModel::set_to_unprepared(void)
{
  _is_prepared = false;
}


inline
void
SemiconductorModel::set_n_dopant(const Dopant& dopant)
{
  _n_dopant = dopant;
  _n_dopant.set_type(Dopant::N_TYPE);
}

inline
void
SemiconductorModel::set_p_dopant(const Dopant& dopant)
{
  _p_dopant = dopant;
  _p_dopant.set_type(Dopant::P_TYPE);
}

inline
void
SemiconductorModel::add_recombination_model(
    DriftDiffusionDefs::RecombinationModel recomb_model)
{
  _recombination |= recomb_model;
}

inline
void
SemiconductorModel::remove_recombination_model(
    DriftDiffusionDefs::RecombinationModel recomb_model)
{
  _recombination &= !recomb_model;
}

inline
void
SemiconductorModel::set_SRH_parameters(double tau_n, double tau_p)
{
  _electron_recombination_time = tau_n;
  _hole_recombination_time = tau_p;
}

inline
double
SemiconductorModel::get_donor_density(void) const
{
  return _n_dopant.get_doping_density();
}

inline
double
SemiconductorModel::get_acceptor_density(void) const
{
  return _p_dopant.get_doping_density();
}


inline
double
SemiconductorModel::get_band_gap(void) const
{
  return _conduction_band.band_edge - _valence_band.band_edge;
}


inline
double
SemiconductorModel::get_DOS_factor(void) const
{
  return _DOS_factor;
}

inline
double
SemiconductorModel::get_conduction_band_edge(void) const
{
  return _conduction_band.band_edge;
}

inline
double
SemiconductorModel::get_valence_band_edge(void) const
{
  return _valence_band.band_edge;
}

inline
void
SemiconductorModel::setup_band_edges(void)
{
  double kT = SimulationOptions::T * Constants::k_B;
  electron_vt = hole_vt = kT;
  
  BandProperties& cb = _conduction_band;
  BandProperties& vb = _valence_band;

  cb.effective_DOS =
    get_DOS_factor() * std::pow(kT * cb.effective_mass, 1.5);

  vb.effective_DOS =
    get_DOS_factor() * std::pow(kT * vb.effective_mass, 1.5);

  conduction_band_edge = cb.band_edge;
  valence_band_edge = vb.band_edge;
}


inline
void
SemiconductorModel::calculate_ionized_donors(double arg_e, double kT,
    double& Nd, double& dNd)
{
  const double arg_max = 100;

  double Ed = _n_dopant.get_ionisation_energy();
  double arg = arg_e + Ed / kT;
  if (arg > arg_max)
  {
    Nd = 0.0;
    dNd = 0.0;
  }
  else
  {
    double tmp = std::exp(arg);
    double g  = _n_dopant.get_g_factor();
    double denom = 1 + g * tmp;

    Nd = _n_dopant.get_doping_density() / denom;
    dNd = -(g * tmp * Nd) / (kT * denom);
  }
}

inline
void
SemiconductorModel::calculate_ionized_acceptors(double arg_h, double kT,
    double& Na, double& dNa)
{
  const double arg_max = 100;

  double Ed = _p_dopant.get_ionisation_energy();
  double arg = arg_h + Ed / kT;
  if (arg > arg_max)
  {
    Na = 0.0;
    dNa = 0.0;
  }
  else
  {
    double tmp = std::exp(arg);
    double g  = _p_dopant.get_g_factor();
    double denom = 1 + g * tmp;

    Na = _p_dopant.get_doping_density() / denom;
    dNa = (g * tmp * Na) / (kT * denom);
  }

}

inline
void
SemiconductorModel::calculate_SRH_recombination(void)
{
  double n  = electron_density;
  double p  = hole_density;
  double dn  = electron_density_derivative;
  double dp  = hole_density_derivative;
  double tn  = _electron_recombination_time;
  double tp  = _hole_recombination_time;
  double ni2 = get_intrinsic_density_squared();
  double ni  = std::sqrt(ni2);
  double denom = tp * (n + ni) + tn * (p + ni);
  //double SRH = (n * p - ni2) / denom;
  double G = ni2 / denom;
  double R = (n * p) / denom;
  double SRH = R - G;
  double a = (p - tp * SRH) * dn / denom;
  double b = (n - tn * SRH) * dp / denom; 
  electron_recombination_rate += SRH;
  electron_recombination_rate_derivatives[1] += a;
  electron_recombination_rate_derivatives[2] += b;
  electron_recombination_rate_derivatives[0] += a + b;
  hole_recombination_rate += SRH;
  hole_recombination_rate_derivatives[1] += a;
  hole_recombination_rate_derivatives[2] += b;
  hole_recombination_rate_derivatives[0] += a + b;
}

inline
void
SemiconductorModel::calculate_direct_recombination(void)
{
  double n  = electron_density;
  double p  = hole_density;
  double dn  = electron_density_derivative;
  double dp  = hole_density_derivative;
  double C  = _direct_rec_param;
  double ni2 = get_intrinsic_density_squared();

  double rec = C * (n * p - ni2);
  double a = C * (dn * p);
  double b = C * (n * dp);
  electron_recombination_rate += rec;
  electron_recombination_rate_derivatives[1] += a;
  electron_recombination_rate_derivatives[2] += b;
  electron_recombination_rate_derivatives[0] += a + b;
  hole_recombination_rate += rec;
  hole_recombination_rate_derivatives[1] += a;
  hole_recombination_rate_derivatives[2] += b;
  hole_recombination_rate_derivatives[0] += a + b;
}
 
inline
void
SemiconductorModel::calculate_exciton_generation(void)
{
  double n  = electron_density;
  double p  = hole_density;
  double dn  = electron_density_derivative;
  double dp  = hole_density_derivative;
  double C  = _exciton_gen_param;
  double ni2 = get_intrinsic_density_squared();

  double rec = C * (n * p - ni2);
  double a = C * (dn * p);
  double b = C * (n * dp);
  electron_recombination_rate += rec;
  electron_recombination_rate_derivatives[1] += a;
  electron_recombination_rate_derivatives[2] += b;
  electron_recombination_rate_derivatives[0] += a + b;
  hole_recombination_rate += rec;
  hole_recombination_rate_derivatives[1] += a;
  hole_recombination_rate_derivatives[2] += b;
  hole_recombination_rate_derivatives[0] += a + b;
}
  
inline
double
SemiconductorModel::get_exciton_generation_rate(void) const
{
  double n  = electron_density;
  double p  = hole_density;
  double ni2 = get_intrinsic_density_squared();

  double rec = _exciton_gen_param * (n * p - ni2);
  return rec;
}
 

inline
void
SemiconductorModel::calculate_Auger_recombination(void)
{
}


template <int coupling>
inline
void
SemiconductorModel::calculate_all(double potential,
    double fermi_e, double fermi_h)
{

  // in this simple model all temperatures are equal
  double kT = electron_vt;
  
  // 1 conduction band
  const BandProperties& cb = _conduction_band;
  const BandProperties& vb = _valence_band;

  double Ec = conduction_band_edge;
  double Ev = valence_band_edge;
  
  // 1.) electron and hole density
  double n = 0, dn = 0, dn2 = 0, dn_over_n = 0, arg_e;
  double p = 0, dp = 0, dp2 = 0, dp_over_p = 0, arg_h;
  if (coupling & DriftDiffusionDefs::ELECTRONS)
  {
    arg_e = (fermi_e + potential - Ec) / kT;
    if (get_statistics() == TiberCad::FERMIDIRAC)
    {
      density_and_derivatives<TiberCad::FERMIDIRAC>(arg_e,
          n, dn, dn2, dn_over_n);
    }
    else
    {
      density_and_derivatives<TiberCad::BOLTZMANN>(arg_e,
          n, dn, dn2, dn_over_n);
    }
  
    double Nc = cb.effective_DOS;
    n *= Nc;
    dn *= Nc / kT;
    dn2 *= Nc / (kT * kT);
    dn_over_n /= kT;

    electron_density = n;
    electron_density_derivative = dn;
  }

  if (coupling & DriftDiffusionDefs::HOLES)
  {
    arg_h = -(fermi_h + potential - Ev) / kT;

    if (get_statistics() == TiberCad::FERMIDIRAC)
    {
      density_and_derivatives<TiberCad::FERMIDIRAC>(arg_h,
          p, dp, dp2, dp_over_p);
    }
    else
    {
      density_and_derivatives<TiberCad::BOLTZMANN>(arg_h,
          p, dp, dp2, dp_over_p);
    }

    double Nv = vb.effective_DOS;
    p *= Nv;
    dp *= -Nv / kT;
    dp2 *= Nv / (kT * kT);
    dp_over_p /= -kT;

    hole_density = p;
    hole_density_derivative = dp;
  }

  // 2.) ionized dopant densities
  double Nd = 0, dNd = 0;
  double Na = 0, dNa = 0;
  if (coupling & DriftDiffusionDefs::ELECTRONS)
  {
    calculate_ionized_donors(arg_e, kT, Nd, dNd);
    ionized_donor_density = Nd;
    ionized_donor_density_derivative = dNd;
  }
  if (coupling & DriftDiffusionDefs::HOLES)
  {
    calculate_ionized_acceptors(arg_h, kT, Na, dNa);
    ionized_acceptor_density = Na;
    ionized_acceptor_density_derivative = dNa;
  }

  // 3.) total charge density
  charge_density = p - n + Nd - Na;
  charge_density_derivatives[0] = dp - dn + dNd - dNa;
  charge_density_derivatives[1] =    - dn + dNd;
  charge_density_derivatives[2] = dp            - dNa;
  
  // 4.) mobilities / conductivities
  // For both statistics:
  // 
  //   mu_n = e * D_n * (1 / n) * (dn / dEf_e)
  //
  // NOTE: kT := kB * T / e includes already e
  if (coupling & DriftDiffusionDefs::ELECTRONS)
  {
    double electron_diffusivity = kT * cb.mobility;
    electron_mobility = electron_diffusivity * dn_over_n;
    electron_conductivity = electron_diffusivity * dn;
    electron_conductivity_derivatives[0] = electron_diffusivity * dn2;
    electron_conductivity_derivatives[1] = electron_diffusivity * dn2;
  }
  if (coupling & DriftDiffusionDefs::HOLES)
  {
    double hole_diffusivity = kT * vb.mobility;
    hole_mobility = -hole_diffusivity * dp_over_p;
    hole_conductivity = -hole_diffusivity * dp;
    hole_conductivity_derivatives[0] = -hole_diffusivity * dp2;
    hole_conductivity_derivatives[2] = -hole_diffusivity * dp2;
  }
  
  electron_recombination_rate = 0;
  electron_recombination_rate_derivatives[0] = 0;
  electron_recombination_rate_derivatives[1] = 0;
  electron_recombination_rate_derivatives[2] = 0;
  hole_recombination_rate = 0;
  hole_recombination_rate_derivatives[0] = 0;
  hole_recombination_rate_derivatives[1] = 0;
  hole_recombination_rate_derivatives[2] = 0;
  // 5.) Recombination
  if (coupling & DriftDiffusionDefs::BOTH)
  {
    if (_recombination & DriftDiffusionDefs::SRH)
      calculate_SRH_recombination();
    if (_recombination & DriftDiffusionDefs::AUGER)
      calculate_Auger_recombination();
    if (_recombination & DriftDiffusionDefs::DIRECT)
      calculate_direct_recombination();

    if (_exciton_gen_param > 1e-56)
      calculate_exciton_generation();
  }
}


#endif //_SEMICONDUCTORMODEL_H_
