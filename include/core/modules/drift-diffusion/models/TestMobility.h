// $Id$

#ifndef _TESTMOBILITY_H_
#define _TESTMOBILITY_H_

#include "MobilityModelInterface.h"


class TBDLEXPORT TestMobility : public MobilityModelInterface
{

  public:

    //! constructor
    TBDLLOCAL TestMobility(void);

    //! Destructor
    virtual ~TestMobility(void);

    //! Create a TestMobility object
    static TestMobility* create(void);

    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void);

    //! \copydoc MobilityModelInterface::get_mobility_derivatives()
    virtual void get_mobility_derivatives(std::vector<double>& dm);

    //! \copydoc MobilityModelInterface::get_derivative_grad_fermi()
    virtual void get_derivative_grad_fermi(RealGradient& dm);


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
    double _mu_max;


    //! The exponent for the temperetaure dependence of _beta
    double _mu_min;


    //! The maximum saturation velocity
    double _E_max;


    //! The driving force to be used
    DrivingForce _force;
      

};

//
// inline methods
// 

inline
TestMobility::TestMobility(void)
  : _mu_max(1100.0),
    _mu_min(100.0),
    _E_max(1e6),
    _force(GRADFERMI)
{
}


inline
TestMobility*
TestMobility::create(void)
{
  return new TestMobility();
}


inline
PhysicalModelInterface*
TestMobility::create_new(void) const
{
  return new TestMobility();
}


inline
void
TestMobility::copy_from(const PhysicalModelInterface* rhs)
{
  MobilityModelInterface::copy_from(rhs);

  const TestMobility* mod =
    dynamic_cast<const TestMobility*>(rhs);
  _mu_max = mod->_mu_max;
  _mu_min = mod->_mu_min;
  _E_max = mod->_E_max;
  _force = mod->_force;
}


inline
TestMobility::~TestMobility(void)
{
}

#endif // _TESTMOBILITY_H_
