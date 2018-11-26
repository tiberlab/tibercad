// $Id$

#ifndef _CONSTANTMOBILITY_H_
#define _CONSTANTMOBILITY_H_

#include "MobilityModelInterface.h"


//! Constant mobility model
/*!
 * The constant mobility model assumes no mobility dependence on
 * doping density or electric field. Only temperature dependence is 
 * included.
 * The mobility is calculated from
 * \f[
 * \mu = \mu_{max} \left(\frac{T}{T_0}\right)^{-\gamma}
 * \f]
 */
class TBDLLOCAL ConstantMobility : public MobilityModelInterface
{

  public:

    //! Destructor
    virtual ~ConstantMobility(void);

    //! Create a ConstantMobility object
    static ConstantMobility* create(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void);

    //! \copydoc MobilityModelInterface::get_mobility_derivatives()
    virtual void get_mobility_derivatives(std::vector<double>& dm);

    //! \copydoc MobilityModelInterface::get_derivative_grad_potential()
    virtual void get_derivative_grad_potential(libMesh::RealGradient& dm);

    //! \copydoc MobilityModelInterface::get_derivative_grad_fermi()
    virtual void get_derivative_grad_fermi(libMesh::RealGradient& dm);


  protected:

    //! constructor
    ConstantMobility(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc MobilityModelInterface::create_new()
    virtual PhysicalModelInterface* create_new(void) const;


  private:

    //! The (constant) mobility
    double mu0_;

    //! The exponent for the temperature dependence
    double exp_;

};

//
// inline methods
// 

inline
ConstantMobility::ConstantMobility(const ModelOptions& options)
  : MobilityModelInterface(options),
    mu0_(1000),
    exp_(1)
{
}


inline
ConstantMobility*
ConstantMobility::create(const ModelOptions& options)
{
  return new ConstantMobility(options);
}


inline
PhysicalModelInterface*
ConstantMobility::create_new(void) const
{
  return new ConstantMobility(get_options());
}


inline
ConstantMobility::~ConstantMobility(void)
{
}

#endif // _CONSTANTMOBILITY_H_
