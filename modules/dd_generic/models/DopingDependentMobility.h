/*  
 * This file is part of the tiberCAD module dd_generic.
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
 * \file DopingDependentMobility.h
 * \brief tiberCAD dd_generic module header.
 *
 * \note This file is part of module dd_generic.
 */


#ifndef TC_DOPINGDEPENDENTMOBILITY_H
#define TC_DOPINGDEPENDENTMOBILITY_H

#include "MobilityModelInterface.h"


//! Doping dependent mobility model
/*!
 * The mobility is assumed to depend on
 * doping density in one of the following ways:
 *
 * \li Formula of Masetti et al. (formula 1):
 * \f[
 * \mu  =  \mu_{min,1}e^{-P_c / N} +
 *  \frac{\mu_{const} - \mu_{min,2}}{1 + (N/C_r)^\alpha} -
 *  \frac{\mu_1}{1 + (C_s/N)^\beta}
 * \f]
 * with \f$\mu_{const}\f$ from the ConstantMobility model.
 *
 * \li Formula of Arora (formula 2):
 * \f{eqnarray*}
 * \mu & = & {\mu_{min}}^\ast +
 * \frac{{\mu_d}^\ast}{1+\left(N/{N_0}^\ast\right)^{\alpha^\ast}} \\
 * \mu_{min}^\ast & = & \mu_{min}\left(\frac{T}{T_0}\right)^{\alpha_m} \\
 * {\mu_d}^\ast & = & \mu_d\left(\frac{T}{T_0}\right)^{\alpha_d} \\
 * {N_0}^\ast & = & N_0\left(\frac{T}{T_0}\right)^{\alpha_N} \\
 * \alpha^\ast & = & \alpha\left(\frac{T}{T_0}\right)^{\alpha_a}
 * \f}
 */
class TC_DLLOCAL DopingDependentMobility : public MobilityModelInterface
{

  public:

    //! Destructor
    virtual ~DopingDependentMobility(void);

    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void) override;

    //! \copydoc MobilityModelInterface::get_derivative_grad_potential()
    virtual void get_derivative_grad_potential(libMesh::RealGradient& dm) override;

    //! \copydoc MobilityModelInterface::get_derivative_grad_fermi()
    virtual void get_derivative_grad_fermi(libMesh::RealGradient& dm) override;

  protected:

    //! constructor
    DopingDependentMobility(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::read_database()
    virtual void read_database(void) override;

    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void) override;

    //! Create submodels
    virtual void prepare_submodels(void) override;


  private:

    //! The formula to be used
    /*!
     * 1 means formula of Masetti et al.
     * 2 means formula of Arora
     */
    int formula_;

    //! The minimum mobility
    /*!
     * is mumin1 for formula 1
     */
    double mumin_;

    //! The temperature coefficient for mumin_
    /*!
     * is mumin2 for formula 1
     */
    double am_;

    //! The difference between maximum and minimum mobility
    /*!
     * is mu1 for formula 1
     */
    double mud_;

    //! The temperature coefficient for mud_
    /*!
     * is Cr for formula 1
     */
    double ad_;

    //! The reference doping density
    /*!
     * is Cs for formula 1
     */
    double N0_;

    //! The temperature coefficient for N0_
    /*!
     * is alpha for formula 1
     */
    double an_;

    //! The exponent
    /*!
     * is beta for formula 1
     */
    double a_;

    //! The temperature coefficient for a_
    /*!
     * is Pc for formula 1
     */
    double aa_;


    //! Constant mobility model
    MobilityModelInterface* const_mob_;

};




//
// inline methods
//

inline
DopingDependentMobility::DopingDependentMobility(const ModelOptions& options)
  : MobilityModelInterface(options),
    formula_(2),
    mumin_(2000),
    am_(-1),
    mud_(6000),
    ad_(-2),
    N0_(7e16),
    an_(4),
    a_(1),
    aa_(0),
    const_mob_(0)
{
}



inline
DopingDependentMobility::~DopingDependentMobility(void)
{
}

#endif // TC_DOPINGDEPENDENTMOBILITY_H
