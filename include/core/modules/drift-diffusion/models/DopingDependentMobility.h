// $Id$

#ifndef _DOPINGDEPENDENTMOBILITY_H_
#define _DOPINGDEPENDENTMOBILITY_H_

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
class DopingDependentMobility : public MobilityModelInterface
{

  public:

    //! Destructor
    virtual ~DopingDependentMobility(void);

    //! Create a DopingDependentMobility object
    static DopingDependentMobility* create(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void);

    //! \copydoc MobilityModelInterface::get_mobility_derivatives()
    virtual void get_mobility_derivatives(std::vector<double>& dm);


  protected:

    //! constructor
    DopingDependentMobility(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void);

    //! Create submodels
    virtual void create_submodels(void);

    //! \copydoc MobilityModelInterface::create_new()
    virtual PhysicalModelInterface* create_new(void) const;


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
DopingDependentMobility*
DopingDependentMobility::create(const ModelOptions& options)
{
  return new DopingDependentMobility(options);
}


inline
PhysicalModelInterface*
DopingDependentMobility::create_new(void) const
{
  return new DopingDependentMobility(get_options());
}


inline
DopingDependentMobility::~DopingDependentMobility(void)
{
  destroy(const_mob_);
}

#endif // _DOPINGDEPENDENTMOBILITY_H_
