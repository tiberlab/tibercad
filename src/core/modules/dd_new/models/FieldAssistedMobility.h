// $Id: FieldAssistedMobility.h 2399 2011-02-23 15:43:57Z maufder $

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
class TBDLLOCAL FieldAssistedMobility : public MobilityModelInterface
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

    //! \copydoc MobilityModelInterface::get_derivative_grad_potential()
    virtual void get_derivative_grad_potential(libMesh::RealGradient& dm);

    //! \copydoc MobilityModelInterface::get_derivative_grad_fermi()
    virtual void get_derivative_grad_fermi(libMesh::RealGradient& dm);


  protected:

    //! constructor
    FieldAssistedMobility(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc MobilityModelInterface::create_new()
    virtual PhysicalModel* create_new(void) const;


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
PhysicalModel*
FieldAssistedMobility::create_new(void) const
{
  return new FieldAssistedMobility(get_options());
}


#endif // _FIELDASSISTEDMOBILITY_H_
