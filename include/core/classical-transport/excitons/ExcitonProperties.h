// $Id$

#ifndef _EXCITONPROPERTIES_H_
#define _EXCITONPROPERTIES_H_

#include "vector_value.h"

#include "SimulationOptions.h"
#include "PhysicalProperties.h"
#include "DriftDiffusionDefs.h"
#include "TiberCad.h"

// GNU scientific library
//#include <gsl/gsl_sf_fermi_dirac.h>

#include <vector>

// forward declarations
class Point;
class Elem;
class DriftDiffusionProperties;

class ExcitonProperties : public PhysicalProperties
{
    
  public:
       
    //! A default (empty) destructor.
    virtual ~ExcitonProperties(void);

    //! Set the statistics to be used
    /*!
     * \param statistics the statistics
     */
    void set_statistics(TiberCad::Statistics statistics);

    //! Get the statistics to be used
    /*!
     * \return the statistics
     */
    TiberCad::Statistics get_statistics(void) const;

    //! (Re-)Initialize for the given element
    /*!
     * \c reinit() calls \c prepare_element_data() which needs to be
     * implemented in derived classes
     *
     * \param elem the current element
     * \param dd_prop a pointer to the semiconductor model
     */
    void reinit(const Elem* elem, DriftDiffusionProperties* dd_prop);
    
    //! The method that will calculate all needed properties
    /*!
     * This method needs to be implemented by a derived class. It has to
     * calculate at least the following set of parameters:
     * 
     * \li the thermal voltage \f$v_T=k_BT/e\f$ for the excitons
     * \li the exciton density and its derivative
     * \li the exciton mobility
     * \li the exciton recombination rate and its derivative
     * \li the exciton generation rate
     * 
     * \param fermi_x the exciton electro-chemical potential
     * \param p the coordinates in real space
     *
     */
    virtual void calculate_all(double fermi_x, const Point& coord) = 0;
      

    //! Get the exciton density
    /*!
     * Get the exciton density as calculated by \c calculate_all(...)
     * 
     * \return the exciton density
     */
    double get_density(void) const
      { return density; };
     
    //! Get the exciton density derivative
    /*!
     * \return the exciton density derivative with respect to the
     * electro-chemical potential
     */
    double get_density_derivative(void) const
      { return density_derivative; };
    

    //! Get the exciton recombination rate
    double get_recombination_rate(void) const
      { return recombination_rate; };
      
    //! Get the exciton recombination rate derivative
    /*!
     * Get \f$\frac{\partial R}{\partial\phi_x}\f$
     */
    double get_recombination_rate_derivative(void) const
        { return recombination_rate_derivative; };
    
    double get_generation_rate(void) const
      { return generation_rate; };
      
    //! Get the exciton mobility
    /*!
     * \return the exciton mobility
     */
    double get_mobility(void) const
      { return mobility; };
      
  protected:
  
    //! The empty constructor.
    ExcitonProperties(void);

    //! This method gets called from reinit()
    /*!
     * It can be used to setup data that is constant in an element, e.g.
     * strain related stuff, band edges.
     */
    virtual void prepare_element_data(void) {};

    //! The element we are currently working on
    const Elem* elem;

    //! The thermal voltage for the electrons
    double exciton_vt;

    //! The density
    double density;

    //! The density derivative
    double density_derivative;

    //! The mobility
    double mobility;

    //! The recombination rate
    double recombination_rate;
    
    //! The derivative of the net electron recombination rate
    double recombination_rate_derivative;
    
    //! The generation rate
    double generation_rate;
    
    //! The band gap
    /*!
     * The band gap is an element data
     */
    double band_gap;

    //! Calculate the density and its derivative
    void calculate_density_and_derivative(double arg, double& density,
        double& derivative) const;
    
    //! Calculate the density for a given argument
    double calculate_density(double arg) const;

    //! Get a reference to the DriftDiffusionProperties
    const DriftDiffusionProperties* get_driftdiffusion_properties(void) const;


  private:

    TiberCad::Statistics _statistics;

    DriftDiffusionProperties* _dd_prop;


};


//
// inline members
//

inline
ExcitonProperties::~ExcitonProperties(void)
{
}

inline
ExcitonProperties::ExcitonProperties(void)
  : PhysicalProperties("ExcitonProperties"),
    elem(NULL),
    _statistics(TiberCad::BOLTZMANN),
    _dd_prop(NULL)
{
}

inline
void
ExcitonProperties::reinit(const Elem* elem, DriftDiffusionProperties* dd_prop)
{
  this->elem = elem;
  this->_dd_prop = dd_prop;
  this->prepare_element_data();
}

inline
void
ExcitonProperties::set_statistics(TiberCad::Statistics statistics)
{
  _statistics = statistics;
}

inline
TiberCad::Statistics
ExcitonProperties::get_statistics(void) const
{
  return _statistics;
}

inline
double
ExcitonProperties::calculate_density(double arg) const
{
  
  const double arg_max = 150;
  const double arg_min = -100;

  double dens;
  if (arg < arg_max)
    dens = std::exp(arg);
  else
    dens = std::exp(arg_max);

  return dens;
}


inline
void
ExcitonProperties::calculate_density_and_derivative(double arg, double& density,
    double& derivative) const
{
  density = calculate_density(arg);
  derivative = density;
}

inline
const DriftDiffusionProperties*
ExcitonProperties::get_driftdiffusion_properties(void) const
{
  return _dd_prop;
}



#endif /* _EXCITONPROPERTIES_H_*/
