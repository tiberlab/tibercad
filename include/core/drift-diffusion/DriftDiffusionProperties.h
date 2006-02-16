// $Id$

#ifndef _DRIFTDIFFUSIONPROPERTIES_H_
#define _DRIFTDIFFUSIONPROPERTIES_H_

//#include "tensor_value.h"
#include "vector_value.h"

#include "PhysicalProperties.h"
#include "TiberCad.h"

#include <vector>

// forward declarations
class Point;
class Elem;

class DriftDiffusionProperties : public PhysicalProperties
{
    
  public:
       
    //! A default (empty) destructor.
    virtual ~DriftDiffusionProperties(void);

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
      

    //! The method that will calculate all needed properties
    /*!
     * This method needs to be implemented by a derived class. It has to
     * calculate at least the following set of parameters:
     * 
     * \li the dielectric tensor
     * \li the total electric polarization
     * \li the electron density and its derivative
     * \li the hole density and its derivative
     * \li the ionized donor density and its derivative
     * \li the ionized acceptor density and its derivative
     * \li the total charge density and its derivatives
     * \li the electron conductivity and its derivatives
     * \li the electron mobility
     * \li the hole mobility
     * \li the hole conductivity and its derivatives
     * \li the net electron recombination rate and its derivatives
     * \li the net hole recombination rate and its derivatives
     * 
     * \param potential the electric potential
     * \param fermi_e the electron electro-chemical potential
     * \param fermi_h the hole electro-chemical potential
     * \param p the coordinates in real space
     * \param elem the pointer to the element which contains \c p
     */
    virtual void calculate_all(double potential,
      double fermi_e, double fermi_h,
      const Point& coord, const Elem* elem) = 0;
      

    //! Get the electron density
    /*!
     * Get the electron density as calculated by \c calculate_all(...)
     * 
     * \return the electron density
     */
    double get_electron_density(void) const
      { return electron_density; };
     
    //! Get the electron density derivative
    /*!
     * \return the electron density derivative with respect to the potential
     */
    //double get_electron_density_derivative(void) const
    //  { return electron_density_derivative; };
    
    //! Get the hole density
    /*!
     * Get the hole density as calculated by \c calculate_all(...)
     * 
     * \return the hole density
     */
    double get_hole_density(void) const
      { return hole_density; };
     
    //! Get the ehole density derivative
    /*!
     * \return the hole density derivative with respect to the potential
     */
    //double get_hole_density_derivative(void) const
    //  { return hole_density_derivative; };
    
    //! Get the ionized donor density
    /*!
     * Get the ionized donor density as calculated by \c calculate_all(...)
     * 
     * \return the ionized donor density
     */
    double get_ionized_donor_density(void) const
      { return ionized_donor_density; };
     
    //! Get the ionized donor density derivative
    /*!
     * \return the ionized donor density derivative with respect to the potential
     */
    //double get_ionized_donor_density_derivative(void) const
    //  { return ionized_donor_density_derivative; };
        
    //! Get the ionized acceptor density
    /*!
     * Get the ionized acceptor density as calculated by \c calculate_all(...)
     * 
     * \return the ionized acceptor density
     */
    double get_ionized_acceptor_density(void) const
      { return ionized_acceptor_density; };
     
    //! Get the ionized acceptor density derivative
    /*!
     * \return the ionized acceptor density derivative with respect to the
     * potential
     */
    //double get_ionized_acceptor_density_derivative(void) const
    //  { return ionized_acceptor_density_derivative; };
    
    
    //! Get the total charge density
    /*!
     * Get the total charge density as calculated by \c calculate_all(...)
     * \f[ \rho = e(p - n + N_D^+ - N_A^-) \f]
     * 
     * \return the total charge density
     */
    double get_charge_density(void) const;
    
    //! Get the derivatives of the charge density
    /*!
     * Get the derivatives of the total charge density with respect to the 
     * electric potential and the two electro-chemical potentials:
     * \f[
     *   \frac{\partial\rho}{\partial\varphi},\;
     *   \frac{\partial\rho}{\partial\phi_n},\;
     *   \frac{\partial\rho}{\partial\phi_p}
     * \f]
     * 
     * \return the derivatives as a const vector reference
     */
    const std::vector<double>& get_charge_density_derivatives(void) const;
    
    //! Get the derivatives of the charge density
    /*!
     * Get the derivatives of the total charge density with respect to the 
     * electric potential one of the two electro-chemical potentials:
     * \f{eqnarray*}
     *   \frac{\partial\rho}{\partial\varphi} & \mathsf{if} i=0, \\
     *   \frac{\partial\rho}{\partial\phi_n} & \mathsf{if} i=1, \\
     *   \frac{\partial\rho}{\partial\phi_p} & \mathsf{if} i=2
     * \f}
     */
    double get_charge_density_derivative(int i) const;
    
    //! Get the net electron recombination rate
    /*!
     * Get \f$R_{net} = R - G\f$ as
     * calculated by \c calculate_all(...)
     */
    double get_net_electron_recombination_rate(void) const
      { return electron_recombination_rate; };
      
    //! Get the net electron recombination rate derivative
    /*!
     * Get \f$frac{\partial R_{net}}{\partial\varphi}\f$
     *
     * \return the derivatives as a const vector reference
     */
    const std::vector<double>&
      get_net_electron_recombination_rate_derivatives(void) const
        { return electron_recombination_rate_derivatives; };
    
    //! Get the net hole recombination rate
    /*!
     * Get \f$R_{net} = R - G\f$ as
     * calculated by \c calculate_all(...)
     *
     */
    double get_net_hole_recombination_rate(void) const
      { return hole_recombination_rate; };
      
    //! Get the net hole recombination rate derivative
    /*!
     * Get \f$frac{\partial R_{net}}{\partial\varphi}\f$
     *
     * \return the derivatives as a const vector reference
     */
    const std::vector<double>&
      get_net_hole_recombination_rate_derivatives(void) const
        { return hole_recombination_rate_derivatives; };

    //! Get the total electric polarization
    /*!
     * The total electric polarization \b P is the sum of the
     * pyroelectric and piezoelectric polarization
     */
    const RealVectorValue& get_total_polarization(void) const
      { return polarization; };

    //! Get the relative permittivity tensor
    //const RealTensorValue& get_relative_permittivity(void) const
    double get_relative_permittivity(void) const
      { return permittivity; };
      
    //! Get the electron conductivity
    /*!
     * \return the electron conductivity \f$\sigma_n = \mu_n n\f$
     */
    double get_electron_conductivity(void) const
      { return electron_conductivity; };
      
    //! Get the hole conductivity
    /*!
     * \return the hole conductivity \f$\sigma_p = \mu_p p\f$
     */
    double get_hole_conductivity(void) const
      { return hole_conductivity; };

    //! Get the electron conductivity derivatives
    const std::vector<double>& get_electron_conductivity_derivatives(void) const
      { return electron_conductivity_derivatives; };
      
    //! Get the hole conductivity derivatives
    const std::vector<double>& get_hole_conductivity_derivatives(void) const
      { return hole_conductivity_derivatives; };

  

  protected:
  
    //! The empty constructor.
    DriftDiffusionProperties(void);

    //! The electron density
    double electron_density;

    //! The electron density derivative
    double electron_density_derivative;

    //! The hole density
    double hole_density;

    //! The hole density derivative
    double hole_density_derivative;
  
    //! The ionized donor density
    double ionized_donor_density;

    //! The ionized donor density derivative
    //double ionized_donor_density_derivative;
        
    //! The ionized acceptor density
    double ionized_acceptor_density;

    //! The ionized acceptor density derivative
    //double ionized_acceptor_density_derivative;

    //! The total charge density
    double charge_density;

    //! The total charge density derivatives
    std::vector<double> charge_density_derivatives;

    //! The electron mobility
    double electron_mobility;

    //! The hole mobility
    double hole_mobility;
    
    //! The electron conductivity
    double electron_conductivity;
    
    //! The derivatives of the electron conductivity
    std::vector<double> electron_conductivity_derivatives;
    
    //! The hole conductivity
    double hole_conductivity;
    
    //! The derivatives of the hole conductivity
    std::vector<double> hole_conductivity_derivatives;
    
    //! The net electron recombination rate
    double electron_recombination_rate;
    
    //! The derivatives of the net electron recombination rate
    std::vector<double> electron_recombination_rate_derivatives;
    
    //! The net hole recombination rate
    double hole_recombination_rate;
    
    //! The derivatives of the net hole recombination rate
    std::vector<double> hole_recombination_rate_derivatives;
    
    //! The total electric polarization
    RealVectorValue polarization;
    
    //! The relative permittivity tensor
    //RealTensorValue permittivity;
    double permittivity;

  private:

    TiberCad::Statistics _statistics;


};


//
// inline members
//

inline
DriftDiffusionProperties::~DriftDiffusionProperties(void)
{
}

inline
DriftDiffusionProperties::DriftDiffusionProperties(void)
  : PhysicalProperties("DriftDiffusionProperties"),
    charge_density_derivatives(3, 0.0),
    electron_conductivity_derivatives(3, 0.0),
    hole_conductivity_derivatives(3, 0.0),
    electron_recombination_rate_derivatives(3, 0.0),
    hole_recombination_rate_derivatives(3, 0.0),
    _statistics(TiberCad::BOLTZMANN)
{
}
    
inline
double
DriftDiffusionProperties::get_charge_density(void) const
{
  return charge_density;
}

inline
const std::vector<double>&
DriftDiffusionProperties::get_charge_density_derivatives(void) const
{
  return charge_density_derivatives;
}

inline
double
DriftDiffusionProperties::get_charge_density_derivative(int i) const
{
  double drho;
  
  switch (i)
  {
    case 1:
      drho = charge_density_derivatives[1];
      break;
    case 2:
      drho = charge_density_derivatives[2];
      break;
    default: // i = 0
      drho = charge_density_derivatives[0];
      break;
  }
  
  return drho;
}

inline
void
DriftDiffusionProperties::set_statistics(TiberCad::Statistics statistics)
{
  _statistics = statistics;
}

inline
TiberCad::Statistics
DriftDiffusionProperties::get_statistics(void) const
{
  return _statistics;
}


#endif /* _DRIFTDIFFUSIONPROPERTIES_H_*/
