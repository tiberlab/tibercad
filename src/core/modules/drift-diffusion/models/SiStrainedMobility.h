// $Id$

#ifndef _SISTRAINEDMOBILITY_H_
#define _SISTRAINEDMOBILITY_H_

#include "MobilityModelInterface.h"
#include "StrainInterface.h"


//! SiStrained mobility model
/*!
 * The constant mobility model assumes no mobility dependence on
 * doping density or electric field. Only temperature dependence is 
 * included.
 * The mobility is calculated from
 * \f[
 * \mu = \mu_{max} \left(\frac{T}{T_0}\right)^{-\gamma}
 * \f]
 */
class TBDLLOCAL SiStrainedMobility : public MobilityModelInterface
{

  public:

    //! Destructor
    virtual ~SiStrainedMobility(void);

    //! Create a SiStrainedMobility object
    static SiStrainedMobility* create(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void);

    //! \copydoc MobilityModelInterface::get_mobility_derivatives()
    virtual void get_mobility_derivatives(std::vector<double>& dm);


  protected:

    //! constructor
    SiStrainedMobility(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc MobilityModelInterface::create_new()
    virtual PhysicalModel* create_new(void) const;


  private:

    //! The strain simulation
    StrainInterface _strain;

    //! The (constant) mobility
    double mu0_;

    //! The exponent for the temperature dependence
    double exp_;

    //! \f$\alpha\f$ for compressive strain
    double _alpha_c;

    //! \f$\beta\f$ for compressive strain
    double _beta_c;

    //! \f$\gamma\f$ for compressive strain
    double _gamma_c;

    //! \f$\alpha\f$ for tensile strain
    double _alpha_t;

    //! \f$\beta\f$ for tensile strain
    double _beta_t;

    //! \f$\gamma\f$ for tensile strain
    double _gamma_t;

};

//
// inline methods
// 




inline
SiStrainedMobility*
SiStrainedMobility::create(const ModelOptions& options)
{
  return new SiStrainedMobility(options);
}


inline
PhysicalModel*
SiStrainedMobility::create_new(void) const
{
  return new SiStrainedMobility(get_options());
}


inline
SiStrainedMobility::~SiStrainedMobility(void)
{
}

#endif // _CONSTANTMOBILITY_H_
