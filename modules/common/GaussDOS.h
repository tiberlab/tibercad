/*  
 * This file is part of the tiberCAD module common.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file GaussDOS.h
 * \brief tiberCAD common module header.
 *
 * \note This file is part of module common.
 */


#ifndef TC_GAUSSDOS_H
#define TC_GAUSSDOS_H


#include "tibercad/physics/misc/DensityOfStates.h"


/*!
* This calculates the G(z,s) function and its derivative with respect to z defined in  G. Paasch, S. Scheinert, Charge carrier density of organics with Gaussian density of states: analytical approximation for the Gauss-Fermi intergral, J. Appl. Phys., 107 (2010) 104501
* \f{equation*}
*	\frac{n}{{N_0 }} = \frac{1}{{\sigma \sqrt {2\pi } }}\int_{ - \infty }^{ + \infty } {dE\exp \left( { - \frac{1}{2}\frac{{\left( {E - E_0 } \right)^2 }}{{\sigma ^2 }}} \right)\frac{1}{{e^{\left( {E - E_F } \right)/kT}  + 1}} \equiv G\left( {\zeta ;s} \right)} 
* \f
* where
* \f{eqnarray*}
* 	\zeta  = \frac{{E_F  - E_0 }}{{kT}} \\
*	s = \frac{\sigma }{{kT}}
* \f
* The G function is calculated according to the following approximation
* If z <= -s^2
* \f{equation*}
* 	G\left( {\zeta ,s} \right) = \exp \left( {\frac{{s^2 }}{2} + \zeta } \right)\frac{1}{{\exp \left[ {K\left( s \right)\left( {\zeta  + s^2 } \right)} \right] + 1}}
* \f
*If z > -s^2
* \f{equation*}
* 	G\left( {\zeta ,s} \right) = \frac{1}{2}\operatorname{erfc} \left( { - \frac{\zeta }{{s\sqrt 2 }}H\left( s \right)} \right)
* \f
*/
class TC_DLEXPORT GaussDOS : public DensityOfStates
{

  public:

    //! Destructor
    virtual ~GaussDOS(void) {};

    //! Creator function
    static GaussDOS* create(const ModelOptions& options);

  protected:

    //! Constructor
    GaussDOS(const ModelOptions& options);

    virtual void read_database(void);

    virtual void do_print_info(void);

    //! Get occupied states and the derivative with respect to phi
    virtual void
    calculate_density_and_derivative(std::vector<double>& result, double E, double Epot, double kT,
                                                                       double kTlattice, const Elem* elem, const Point& p) const;
   
    virtual void do_init(void);

  private:

    //Gauss parameters
    double _E0;
    double _sigma;

    //Total density parameter
    double _N0;

    /*!
    * This calculates the K(s) function defined in  G. Paasch, S. Scheinert, Charge carrier density of organics with Gaussian density of states: analytical approximation for the Gauss-Fermi intergral, J. Appl. Phys., 107 (2010) 104501
    * \f{equation*}
    * 	K\left( x \right) = 2\left\{ {1 - \frac{H}{x}\sqrt {\frac{2}{\pi }} \exp \left[ {\frac{1}{2}x^2 \left( {1 - H^2 } \right)} \right]} \right\}
    * \f
    */
    double K (double x, double h) const;

    /*!
    * This calculates the H(s) function defined in  G. Paasch, S. Scheinert, Charge carrier density of organics with Gaussian density of states: analytical approximation for the Gauss-Fermi intergral, J. Appl. Phys., 107 (2010) 104501
    * \f{equation*}
    * 	H\left( x \right) = \frac{{\sqrt 2 }}{x}erfc^{ - 1} \left[ {\exp \left( { - \frac{{x^2 }}{2}} \right)} \right]
    * \f
    */
    double H (double x) const;

    /*! The erfc = 1 - erf is calculated using a Taylor series expansion for x <= 4 (corresponding to about 5.66 standard deviations) 
    * and for x > 4 using the asymptotic expansion 7.1.26 given in Abramowitz and Stegun (http://people.math.sfu.ca/~cbm/aands/page_299.htm)
    */
    double erfc (double x) const;	

    /*!
    * The inverse erfc is calculated for  x > 0.002 as 
    *\f{eqnarray*}
    * 	erfc^{ - 1} \left( x \right) =  - \sum\limits_{n = 0}^{ + \infty } {\frac{{c_n }}{{2n + 1}}\left[ {\frac{{\sqrt \pi  }}{2}\left( {x - 1} \right)}\right]^{2n + 1} }  \\
    *	c_0 = 1 \\
    *	c_n  = \sum\limits_{m = 0}^{n - 1} {\frac{{c_m c_{n - 1 - m} }}{{\left( {m + 1} \right)\left( {2m + 1} \right)}}} 
    *\f
    * (see http://functions.wolfram.com/GammaBetaErf/InverseErfc/introductions/ProbabilityIntegrals/ShowAll.html)
    * and for x <= 0.002 as the following asymptotic expansion ( JR Philip, The Function Inverfc, Australian Journal of Physics 13 (1960) 13 )
    *\f{equation*}
    * 	erfc^{ - 1} \left( x \right) \sim \sqrt { - \log \left( x \right) - \frac{1}{2}\log \left( { - \pi \log \left( x \right) - \frac{1}{2}\log \left( { - \pi \log \left( x \right)} \right)} \right)} 
    *\f
    */
    double inverfc (double x) const;
};

//
// inline methods
//

inline
GaussDOS*
GaussDOS::create(const ModelOptions& options)
{
  return new GaussDOS(options);
}


#endif // TC_GAUSSDOS_H
