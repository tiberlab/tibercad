// $Id$

#ifndef _SEMICONDUCTORMODEL_H_
#define _SEMICONDUCTORMODEL_H_

#include "SimulationOptions.h"
#include "Dopant.h"
#include "DriftDiffusionProperties.h"

#include <vector>
#include <string>
#include <cmath>


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

    void set_mobilities(double mu_e, double mu_h)
    {
      _e_mobility = mu_e;
      _h_mobility = mu_h;
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
      double fermi_e, double fermi_h, const Point& coord);
    
    //! Get the exciton generation rate
    double get_exciton_generation_rate(void) const;


    void print_info(void) const;

    virtual void get_net_recombination_rates(std::vector<double>& rates);
    
  protected:

    //! Get the physical semiconductor model
    /*!
     * Derived classes will need to access the physical model, e.g. to
     * set the strain.
     */
    DDsemiconductor* get_physical_model(void)
      { return _bulk_model; };

    //! Extract the band properties from _bulk_model
    /*!
     * This method looks for the band extrema and puts the effective
     * mass, band edges etc. into the BandProperties structure
     */
    void extract_band_properties(void);

    //! \copydoc DriftDiffusionProperties::prepare_element_data()
    virtual void prepare_element_data(void);

    //! Calculate Shockley-Read-Hall recombination
    void calculate_SRH_recombination(void);

    //! Calculate Shockley-Read-Hall recombination
    double get_SRH_recombination(void);

    //! Calculate Auger recombination
    void calculate_Auger_recombination(void);

    //! Calculate Auger recombination
    double get_Auger_recombination(void);

    //! Calculate direct recombination
    void calculate_direct_recombination(void);

    //! Calculate direct recombination
    double get_direct_recombination(void);

    //! Calculate exciton generation
    void calculate_exciton_generation(void);

    //! Set the object to unprepared state
    void set_to_unprepared(void);

  private:

    typedef DriftDiffusionProperties Parent;

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

    double _e_mobility;
    double _h_mobility;


    std::string _filename;

    
    //! The recombination models used
    int _recombination;

    double _electron_recombination_time;
    double _hole_recombination_time;

    double _direct_rec_param;

    double _exciton_gen_param;

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
double
SemiconductorModel::get_SRH_recombination(void)
{
  double n  = electron_density;
  double p  = hole_density;
  double tn  = _electron_recombination_time;
  double tp  = _hole_recombination_time;
  double ni2 = get_intrinsic_density_squared();
  double ni  = std::sqrt(ni2);
  double denom = tp * (n + ni) + tn * (p + ni);
  double G = ni2 / denom;
  double R = (n * p) / denom;
  return (R - G);
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
double
SemiconductorModel::get_direct_recombination(void)
{
  double n  = electron_density;
  double p  = hole_density;
  double C  = _direct_rec_param;
  double ni2 = get_intrinsic_density_squared();

  double rec = C * (n * p - ni2);
  return rec;
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

inline
double
SemiconductorModel::get_Auger_recombination(void)
{
  return 0;
}


inline
void
SemiconductorModel::get_net_recombination_rates(std::vector<double>& rates)
{
  rates.resize(3);

  if (_recombination & DriftDiffusionDefs::SRH)
    rates[0] = get_SRH_recombination();
  if (_recombination & DriftDiffusionDefs::DIRECT)
    rates[1] = get_direct_recombination();
  if (_recombination & DriftDiffusionDefs::AUGER)
    rates[2] = get_Auger_recombination();
}


#endif //_SEMICONDUCTORMODEL_H_
