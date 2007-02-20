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
 * \mu  =  \mu_{min,1}\exp(-P_c / N) +
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
    static DopingDependentMobility* create(void);

    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void);

    //! \copydoc MobilityModelInterface::get_mobility_derivatives()
    virtual void get_mobility_derivatives(std::vector<double>& dm);


  protected:

    //! constructor
    DopingDependentMobility(void);

    //! \copydoc MobilityModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc MobilityModelInterface::create_new()
    virtual PhysicalModelInterface* create_new(void) const;

    //! \copydoc MobilityModelInterface::copy_from()
    virtual void copy_from(const PhysicalModelInterface* rhs);

    /*! \copydoc MobilityModelInterface::calculate_VCA() */
    virtual void calculate_VCA(const PhysicalModelInterface* comp_A,
        const PhysicalModelInterface* comp_B, double xa);

  private:

    //! The formula to be used
    /*!
     * 1 means formula of Masetti et al.
     * 2 means formula of Arora
     */
    int _formula;

    //! The minimum mobility
    /*!
     * is mumin1 for formula 1
     */
    double _mumin;

    //! The temperature coefficient for _mumin
    /*!
     * is mumin2 for formula 1
     */
    double _am;

    //! The difference between maximum and minimum mobility
    /*!
     * is mu1 for formula 1
     */
    double _mud;

    //! The temperature coefficient for _mud
    /*!
     * is Cr for formula 1
     */
    double _ad;

    //! The reference doping density
    /*!
     * is Cs for formula 1
     */
    double _N0;

    //! The temperature coefficient for _N0
    /*!
     * is alpha for formula 1
     */
    double _an;

    //! The exponent
    /*!
     * is beta for formula 1
     */
    double _a;

    //! The temperature coefficient for _a
    /*!
     * is Pc for formula 1
     */
    double _aa;


    //! Constant mobility model
    MobilityModelInterface* _const_mob;

};




//
// inline methods
// 

inline
DopingDependentMobility::DopingDependentMobility(void)
  : _formula(2),
    _mumin(2000),
    _am(-1),
    _mud(6000),
    _ad(-2),
    _N0(7e16),
    _an(4),
    _a(1),
    _aa(0)
{
}


inline
DopingDependentMobility*
DopingDependentMobility::create(void)
{
  return new DopingDependentMobility();
}


inline
PhysicalModelInterface*
DopingDependentMobility::create_new(void) const
{
  return new DopingDependentMobility();
}


inline
DopingDependentMobility::~DopingDependentMobility(void)
{
}

#endif // _DOPINGDEPENDENTMOBILITY_H_
