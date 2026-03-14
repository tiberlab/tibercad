/*  
 * This file is part of the tiberCAD module driftdiffusion.
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
 * \file HoppingMobility.h
 * \brief tiberCAD driftdiffusion module header.
 *
 * \note This file is part of module driftdiffusion.
 */



#ifndef TC_HOPPINGMOBILITY_H
#define TC_HOPPINGMOBILITY_H

#include "MobilityModelInterface.h"

/*! See Pasveer et al. Unified Description of Charge-Carrier Mobilities in Disordered Semiconducting Polymers, Phys. Rev. Lett. 94 (2005) 206601
* referring to G. Paasch, S. Scheinert, Charge carrier density of organics with Gaussian density of states: analytical approximation for the Gauss-Fermi integral, J. Appl. Phys., 107 (2010) 104501
* Ef: Fermi level;
*E0, sigma: gaussian parameters;
*N0: density of sites;
*F: electric field;
*e: electron charge;
* \f{eqnarray*}
*	\mu = \mu \left(s , G \right) f \left( s, F \right) \\
*	\mu \left(s , G \right) = \mu_0 \left( s \right) \exp \left\[ \frac{1}{2} \left( s^2 - s \right) \left(2G \right)^\delta \right\] \\
*	\delta = 2 \frac{ \ln \left( s^2 - s \right) - \ln \left( \ln 4 \right)}{s^2} \\
*	\mu_0 \left( s \right) = \mu_0 c_1 \exp \left( - c_2 s^2 \right) \\
*	G = G \left(\zeta, s \right) \\
*	\zeta = \frac{E_f - E_0}{kT} \\
*	s = \frac{\sigma}{kT} \\
*	c_1 = 1.8 \cdot 10^{-9} \\
*	c_2 = 0.42 \\
*	\mu_0 = \left(\frac{1}{N_0} \right)^{\frac{2}{3}} \frac{\nu_0 e}{\sigma} \\
*	f \left( s, F \right) = \exp \left\{ 0.44 \left( s^{\frac{3}{2}} - 2.2 \right) \left\[ \sqrt{1 + 0.8 \left(F \frac{e}{\sigma N_0^{\frac{1}{3}}} \right)^2 } - 1 \right\] \right\} 
* \f
*/
class TC_DLLOCAL HoppingMobility : public MobilityModelInterface
{

  public:

    //! Destructor
    virtual ~HoppingMobility(void);

    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void);

    //! \copydoc MobilityModelInterface::get_mobility_derivatives()
    virtual void get_mobility_derivatives(std::vector<double>& dm);

    //! \copydoc MobilityModelInterface::get_derivative_potential()
    virtual double get_derivative_potential(void);

    //! \copydoc MobilityModelInterface::get_derivative_grad_potential()
    virtual void get_derivative_grad_potential(libMesh::RealGradient& dm);

    //! \copydoc MobilityModelInterface::get_derivative_grad_fermi()
    virtual void get_derivative_grad_fermi(libMesh::RealGradient& dm);

  protected:

    //! constructor
    HoppingMobility(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void);

  private:

    //! Attempt to jump frequency according to the Miller Abrahams model
    double _nu0;

    //! Site density
    double _N0;

    //! Variance of gaussian distribution of states
    double _sigma;

};

//
// inline methods
//

inline
HoppingMobility::HoppingMobility(const ModelOptions& options)
  : MobilityModelInterface(options),
    _nu0(1e13),
    _N0(1e21),
    _sigma(0.1)
{
}


inline
HoppingMobility::~HoppingMobility(void)
{
}

#endif // TC_HOPPINGMOBILITY_H
