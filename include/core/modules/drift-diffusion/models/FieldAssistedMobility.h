// $Id$

#ifndef _FIELDASSISTEDMOBILITY_H_
#define _FIELDASSISTEDMOBILITY_H_

#include "MobilityModelInterface.h"


//! Field-assisted mobility model for unordered systems
/*!
 * The mobility is calculated as
 * \f[
 * \mu = \mu_0 e^{\sqrt{|E|/E_0}}
 * \f]
 */
class FieldAssistedMobility : public MobilityModelInterface
{

  public:

    //! Destructor
    virtual ~FieldAssistedMobility(void) {};

    //! Create a FieldAssistedMobility object
    static FieldAssistedMobility* create(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void);

    //! \copydoc MobilityModelInterface::get_mobility_derivatives()
    virtual void get_mobility_derivatives(std::vector<double>& dm);

    //! \copydoc MobilityModelInterface::get_derivative_grad_fermi()
    virtual void get_derivative_grad_fermi(RealGradient& dm);


  protected:

    //! constructor
    FieldAssistedMobility(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc MobilityModelInterface::create_new()
    virtual PhysicalModelInterface* create_new(void) const;


  private:

    enum DrivingForce
    {
      EFIELD,
      GRADFERMI,
    };


    //! The zero-field mobility
    double _mu0;


    //! The critical field strength
    double _E0;


    //! The driving force to be used
    DrivingForce _force;


};

//
// inline methods
//

inline
FieldAssistedMobility::FieldAssistedMobility(const ModelOptions& options)
  : MobilityModelInterface(options),
    _mu0(0.0054),
    _E0(3e5),
    _force(EFIELD)
{
}


inline
FieldAssistedMobility*
FieldAssistedMobility::create(const ModelOptions& options)
{
  return new FieldAssistedMobility(options);
}


inline
PhysicalModelInterface*
FieldAssistedMobility::create_new(void) const
{
  return new FieldAssistedMobility(get_options());
}


#endif // _FIELDASSISTEDMOBILITY_H_
