// $Id$

#ifndef _BODYFORCEMODEL_H_
#define _BODYFORCEMODEL_H_

#include "PhysicalModelInterface.h"

#include "tensor_value.h"
#include "vector_value.h"

class Elem;
class Point;

using namespace std;

//! The base class for body force models
class BodyForceModel : public PhysicalModelInterface
{

  public:

    //! Destructor
    virtual ~BodyForceModel(void) {};

    //! Creator function
    static BodyForceModel* create(const ModelOptions& options);

    const RealGradient& get_force_source(void) const;

    const RealTensor& get_stress_source(void) const;

    const RealTensor& get_strain_source(void) const;


    //! Calculate local body force
    /*!
     * \param elem pointer to the current element
     * \param point the coordinates in the reference element
     */
    virtual void calculate(const Elem* elem, const Point& point){};

  protected:

    //! Constructor
    BodyForceModel(const ModelOptions& options);

    void set_force_source(const RealGradient& force_source);

    void set_strain_source(const RealTensor& strain_source);

    void set_stress_source(const RealTensor& stress_source);



  private:

    RealGradient _force_source;

    RealTensor _strain_source;

    RealTensor _stress_source;

};

inline
const RealGradient&
BodyForceModel::get_force_source(void) const
{
  return _force_source;
}

inline
const RealTensor&
BodyForceModel::get_strain_source(void) const
{
  return _strain_source;
}

inline
const RealTensor&
BodyForceModel::get_stress_source(void) const
{
  return _stress_source;
}

inline 
void 
BodyForceModel::set_force_source(const RealGradient& force_source)
{
  _force_source = force_source;
}

inline 
void 
BodyForceModel::set_strain_source(const RealTensor& strain_source)
{
  _strain_source = strain_source;
}


inline 
void 
BodyForceModel::set_stress_source(const RealTensor& stress_source)
{
  _stress_source = stress_source;
}


inline
BodyForceModel::BodyForceModel(const ModelOptions& options) :
  PhysicalModelInterface(options),
  _force_source(0),
  _strain_source(0),
  _stress_source(0)
{
}

#endif // _THERMALCONDUCTIVITYMODEL_H_
