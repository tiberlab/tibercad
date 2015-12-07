// $Id$

#ifndef _FIELDDEPENDENTMOBILITY_H_
#define _FIELDDEPENDENTMOBILITY_H_

#include "MobilityModelInterface.h"


//! FieldDependent mobility model
/*!
 * The constant mobility model assumes no mobility dependence on
 * doping density or electric field. Only temperature dependence is
 * included.
 * The mobility is calculated from
 * \f[
 * \mu = \mu_{max} \left(\frac{T}{T_0}\right)^{-\gamma}
 * \f]
 */
class TBDLLOCAL FieldDependentMobility : public MobilityModelInterface
{

  public:

    //! Destructor
    virtual ~FieldDependentMobility(void);

    //! Create a FieldDependentMobility object
    static FieldDependentMobility* create(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void);

    //! \copydoc MobilityModelInterface::get_derivative_grad_potential()
    virtual void get_derivative_grad_potential(RealGradient& dm);

    //! \copydoc MobilityModelInterface::get_derivative_grad_fermi()
    virtual void get_derivative_grad_fermi(RealGradient& dm);


  protected:

    //! constructor
    FieldDependentMobility(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void);

    //! Create low field mobility model
    virtual void prepare_submodels(void);



  private:

    enum DrivingForce
    {
      EFIELD,
      GRADFERMI,
      FIELDPARAM
    };


    //! The exponent
    double _beta;


    //! The exponent for the temperetaure dependence of _beta
    double _betaexp;


    //! The maximum saturation velocity
    double _vsat0;


    //! For the temparature dependence of v_sat
    double _vsat_b;


    //! The minimum of vsat for formula 2
    double _vsat_min;


    //! The formula to be used
    int _vsat_formula;


    //! The low-field mobility
    MobilityModelInterface* _low_field_mob;


    //! The driving force to be used
    DrivingForce _force;


    //! A damping parameter
    double _damping;


};

//
// inline methods
//

inline
FieldDependentMobility::FieldDependentMobility(const ModelOptions& options)
  : MobilityModelInterface(options),
    _beta(1),
    _betaexp(0.0),
    _vsat0(1.13e7),
    _vsat_b(1),
    _vsat_min(5e5),
    _vsat_formula(1),
    _low_field_mob(NULL),
    _force(EFIELD),
    _damping(1e9)
{
}


inline
FieldDependentMobility*
FieldDependentMobility::create(const ModelOptions& options)
{
  return new FieldDependentMobility(options);
}



inline
FieldDependentMobility::~FieldDependentMobility(void)
{
}

#endif // _FIELDDEPENDENTMOBILITY_H_
