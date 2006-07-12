// $Id$

#ifndef _SIMPLESEMICONDUCTORMODEL_H_
#define _SIMPLESEMICONDUCTORMODEL_H_

#include "SimulationOptions.h"
#include "Dopant.h"
#include "DriftDiffusionProperties.h"
#include "Constants.h"

#include <vector>
#include <cmath>

// forward declarations
//class Point;
class Elem;

//! A simple semiconductor model
/*!
 * This simple model uses manually provided band properties
 */
class SimpleSemiconductorModel : public DriftDiffusionProperties
{
  public:

    //! The default constructor
    SimpleSemiconductorModel(void);

    //! The destructor
    virtual ~SimpleSemiconductorModel(void) {};

    //! Set conduction band properties
    /*!
     * The DOS effective mass has to include the spin degeneracy.
     */
    void set_conduction_band_properties(double band_edge,
        double effective_mass, double mobility);

    //! Set valence band properties
    /*!
     * The DOS effective mass has to include the spin degeneracy.
     */
    void set_valence_band_properties(double band_edge,
        double effective_mass, double mobility);
    
    //! Set the relative permittivity
    void set_relative_permittivity(double epsilon_r);
    
    virtual void read_database(const Dummy&) {};

    /*! \copydoc DriftDiffusionProperties::calculate_all()
     * 
     * This implementation models the most simple semiconductor equations
     */
    virtual void calculate_all(double potential,
      double fermi_e, double fermi_h, const Point& coord);
    
  protected:

    //! \copydoc DriftDiffusionProperties::prepare_element_data()
    /*!
     * This implementation calculates the effective density of states
     */
    virtual void prepare_element_data(void);

  private:

    typedef DriftDiffusionProperties Parent;
    
    SimpleSemiconductorModel(const SimpleSemiconductorModel& model);
    
    //! \c true if equilibrium properties are calculated
    bool _is_prepared;

};


//
// inline member functions
//


inline
void
SimpleSemiconductorModel::set_conduction_band_properties(double band_edge,
    double effective_mass, double mobility)
{
  get_conduction_band().band_edge = band_edge;
  get_conduction_band().effective_mass = effective_mass;
  electron_mobility = mobility;
}

inline
void
SimpleSemiconductorModel::set_valence_band_properties(double band_edge,
    double effective_mass, double mobility)
{
  get_valence_band().band_edge = band_edge;
  get_valence_band().effective_mass = effective_mass;
  hole_mobility = mobility;
}

inline
void
SimpleSemiconductorModel::set_relative_permittivity(double epsilon_r)
{
  permittivity = epsilon_r;
}


#endif //_SIMPLESEMICONDUCTORMODEL_H_
