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
    SimpleSemiconductorModel(void);

    //! The destructor
    virtual ~SimpleSemiconductorModel(void) {};

    //! This method creates a SimpleSemiconductorModel object
    static SimpleSemiconductorModel* create(void);

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

    //! \copydoc DriftDiffusionProperties::copy_from()
    virtual void copy_from(const PhysicalModelInterface* rhs);


    //! \copydoc DriftDiffusionProperties::prepare_element_data()
    /*!
     * This implementation calculates the effective density of states
     */
    virtual void prepare_element_data(void);

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
SimpleSemiconductorModel::create(void)
{
  return new SimpleSemiconductorModel();
}


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


inline
PhysicalModelInterface*
SimpleSemiconductorModel::create_new(void) const
{
  return new SimpleSemiconductorModel();
}


inline
void
SimpleSemiconductorModel::copy_from(const PhysicalModelInterface* rhs)
{
  Parent::copy_from(rhs);

  const SimpleSemiconductorModel* mod =
    dynamic_cast<const SimpleSemiconductorModel*>(rhs);
  
  is_prepared_ = mod->is_prepared_;
  
  get_conduction_band().band_edge = mod->get_conduction_band().band_edge;
  get_conduction_band().effective_mass =
    mod->get_conduction_band().effective_mass;
  get_valence_band().band_edge = mod->get_valence_band().band_edge;
  get_valence_band().effective_mass = mod->get_valence_band().effective_mass;
}



#endif //_SIMPLESEMICONDUCTORMODEL_H_
