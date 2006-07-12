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
 * This model calculates band properties using the \c DDsemiconductor
 * interface. It is not intended for use in strained structures.
 * 
 */
class SemiconductorModel : public DriftDiffusionProperties
{
  public:

    //! The default constructor
    SemiconductorModel(void);

    //! The destructor
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

    //! Set exciton generation rate parameter
    void set_exciton_generation_rate_parameter(double g)
      { _exciton_gen_param = g; };

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

    //! Calculate exciton generation
    void calculate_exciton_generation(void);

    //! Set the object to unprepared state
    void set_to_unprepared(void);

  private:

    typedef DriftDiffusionProperties Parent;
    
    SemiconductorModel(const SemiconductorModel& model);

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
 
#endif //_SEMICONDUCTORMODEL_H_
