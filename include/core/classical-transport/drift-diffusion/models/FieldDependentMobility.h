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
class FieldDependentMobility : public MobilityModelInterface
{

  public:

    //! constructor
    FieldDependentMobility(void);

    //! Destructor
    virtual ~FieldDependentMobility(void);

    //! Create a FieldDependentMobility object
    static FieldDependentMobility* create(void);

    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void);

    //! \copydoc MobilityModelInterface::get_mobility_derivatives()
    virtual void get_mobility_derivatives(std::vector<double>& dm);


  protected:

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

    enum DrivingForce
    {
      EFIELD,
      GRADFERMI,
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
      

};

//
// inline methods
// 

inline
FieldDependentMobility::FieldDependentMobility(void)
  : _beta(1),
    _betaexp(0.0),
    _vsat0(1.13e7),
    _vsat_b(1),
    _vsat_min(5e5),
    _vsat_formula(1),
    _low_field_mob(NULL),
    _force(EFIELD)
{
}


inline
FieldDependentMobility*
FieldDependentMobility::create(void)
{
  return new FieldDependentMobility();
}


inline
PhysicalModelInterface*
FieldDependentMobility::create_new(void) const
{
  return new FieldDependentMobility();
}


inline
void
FieldDependentMobility::copy_from(const PhysicalModelInterface* rhs)
{
  MobilityModelInterface::copy_from(rhs);

  const FieldDependentMobility* mod =
    dynamic_cast<const FieldDependentMobility*>(rhs);
  _beta = mod->_beta;
  _vsat0 = mod->_vsat0;
  _vsat_b = mod->_vsat_b;
  _vsat_formula = mod->_vsat_formula;
}


inline
FieldDependentMobility::~FieldDependentMobility(void)
{
  PhysicalModelInterface::destroy(_low_field_mob);
}

#endif // _FIELDDEPENDENTMOBILITY_H_
