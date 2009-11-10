// $Id$

#ifndef _PYROPOLARIZATION_H_
#define _PYROPOLARIZATION_H_

#include "PhysicalModelInterface.h"
#include "SimulationOptions.h"
#include "TypeDefs.h"
#include "tensor.h"

class Elem;
class Point;

//! A base class for pyropolarization
class PyroPolarization : public PhysicalModelInterface
{

  public:

    //! Get the pyropolarization as Tensor1
    const Tensor1& get_polarization(void) const;

    //! Create for material mat
    static PyroPolarization* create(const Material* mat);

    //! The create method
    static PyroPolarization* create(void);

    //! Calculate the polarization
    /*!
     * Pyropolarization models will most probably depend on temperature, so we pass
     * it (defaulting to the simulation temperature). The element and point are
     * passed such that the model could implement any type of dependence.
     *
     * Calls do_calculate_polarization()
     */
    void calculate_polarization(const Elem* elem, const Point& p,
        double temperature = SimulationOptions::temperature);


  protected:

    //! Constructor
    PyroPolarization(void);

    //! Destructor
    virtual ~PyroPolarization(void);

    
    //! Get a writable reference to the polarization
    Tensor1& get_polarization_vector(void);
    
    //! Set the polarization in the calculation system
    void set_polarization(const Tensor1& polarization);

    //! Calculates the polarization
    /*!
     * This is called from calculate_polarization() and can be used
     * to implement e.g. temperature dependent models
     */
    virtual void do_calculate_polarization(const Elem* elem, const Point& p,
        double temperature);

  
    //! Create new PyroPolarization
    virtual PhysicalModelInterface* create_new(void) const;
  

  private:

    //! The pyropolarization in the calculation system
    Tensor1 _polarization;

};


//
// inline methods
// 


inline
PyroPolarization::PyroPolarization(void)
  : _polarization(0)
{
  set_name("PyroPolarization");
}


inline
PyroPolarization::~PyroPolarization(void)
{
}


inline
const Tensor1&
PyroPolarization::get_polarization(void) const
{
  return _polarization;
}

inline
void
PyroPolarization::set_polarization(const Tensor1& polarization)
{
  _polarization = polarization;
}

inline
void
PyroPolarization::calculate_polarization(const Elem* elem, const Point& p,
        double temperature)
{
  do_calculate_polarization(elem, p, temperature);
}

inline
void
PyroPolarization::do_calculate_polarization(const Elem* elem, const Point& p,
        double temperature)
{
  ignore_unused_variable(elem);
  ignore_unused_variable(&p);
  ignore_unused_variable(temperature);
}

#endif // _PYROPOLARIZATION_H_
