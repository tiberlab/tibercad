// $Id: InterpolatedMobility.h 4135 2015-09-25 10:19:38Z maufder $

#ifndef _INTERPOLATEDMOBILITY_H_
#define _INTERPOLATEDMOBILITY_H_

#include "MobilityModelInterface.h"

class SimulationInterface;

//! Interpolated mobility model

class TBDLLOCAL InterpolatedMobility : public MobilityModelInterface
{

  public:

    //! Destructor
    virtual ~InterpolatedMobility(void);

    //! Create a InterpolatedMobility object
    static InterpolatedMobility* create(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::get_mobility()
    virtual double get_mobility(void);

    //! \copydoc MobilityModelInterface::get_derivative_potential()
    virtual double get_derivative_potential(void);

    //! \copydoc MobilityModelInterface::get_mobility_derivatives()
    virtual void get_mobility_derivatives(std::vector<double>& dm);

    //! \copydoc MobilityModelInterface::get_derivative_grad_potential()
    virtual void get_derivative_grad_potential(libMesh::RealGradient& dm);

    //! \copydoc MobilityModelInterface::get_derivative_grad_fermi()
    virtual void get_derivative_grad_fermi(libMesh::RealGradient& dm);


  protected:

    //! constructor
    InterpolatedMobility(const ModelOptions& options);

    //! \copydoc MobilityModelInterface::read_database()
    virtual void read_database(void);

    //! \copydoc MobilityModelInterface::do_init()
    virtual void do_init(void);

    //! \copydoc MobilityModelInterface::create_new()
    virtual PhysicalModel* create_new(void) const;


  private:

  //! The interpolation simulation to use
  SimulationInterface* _interpolation_sim;

  //! Model name
  std::string _model_name;

  //Variables
  std::string _temperature_var;
  std::string _field_var;
  std::string _density_var;

  ID _model_id;
  ID _temperature_id;
  ID _field_id;
  ID _density_id;
  
};

//
// inline methods
// 

inline
InterpolatedMobility::InterpolatedMobility(const ModelOptions& options)
  : MobilityModelInterface(options)
{
}


inline
InterpolatedMobility*
InterpolatedMobility::create(const ModelOptions& options)
{
  return new InterpolatedMobility(options);
}


inline
PhysicalModel*
InterpolatedMobility::create_new(void) const
{
  return new InterpolatedMobility(get_options());
}


inline
InterpolatedMobility::~InterpolatedMobility(void)
{
}

#endif // _INTERPOLATEDMOBILITY_H_
