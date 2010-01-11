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
    /*!
     * A derived class can call this constructor to set a different
     * model name
     */
    SimpleSemiconductorModel(const ModelOptions& options);

    //! The destructor
    virtual ~SimpleSemiconductorModel(void) {};

    //! This method creates a SimpleSemiconductorModel object
    static SimpleSemiconductorModel* create(const ModelOptions& options);

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
    

  protected:

    //! \copydoc DriftDiffusionProperties::do_init()
    virtual void do_init();

    //! \copydoc DriftDiffusionProperties::create_new()
    virtual PhysicalModelInterface* create_new(void) const;


    //! \copydoc DriftDiffusionProperties::prepare_element_data()
    /*!
     * This implementation calculates the effective density of states
     */
    virtual void prepare_element_data(void);

    /*! \copydoc PhysicalModelInterface::do_print_info() */
    virtual void do_print_info(void);


  private:

    typedef DriftDiffusionProperties Parent;

    SimpleSemiconductorModel(const SimpleSemiconductorModel& model);
    SimpleSemiconductorModel& operator=(const SimpleSemiconductorModel& model);
    
    //! \c true if equilibrium properties are calculated
    bool is_prepared_;

};


//
// inline member functions
//

inline
SimpleSemiconductorModel*
SimpleSemiconductorModel::create(const ModelOptions& options)
{
  return new SimpleSemiconductorModel(options);
}


inline
void
SimpleSemiconductorModel::set_conduction_band_properties(double band_edge,
    double effective_mass, double mobility)
{
  get_conduction_band().band_edge = band_edge;
  get_conduction_band().effective_mass = effective_mass;
  get_pd().electron_mobility = mobility;
}


inline
void
SimpleSemiconductorModel::set_valence_band_properties(double band_edge,
    double effective_mass, double mobility)
{
  get_valence_band().band_edge = band_edge;
  get_valence_band().effective_mass = effective_mass;
  get_pd().hole_mobility = mobility;
}


inline
void
SimpleSemiconductorModel::set_relative_permittivity(double epsilon_r)
{
  permittivity = epsilon_r;
}


inline
PhysicalModelInterface*
SimpleSemiconductorModel::create_new(void) const
{
  return new SimpleSemiconductorModel(get_options());
}



#endif //_SIMPLESEMICONDUCTORMODEL_H_
