// $Id: PoissonBoundaryModel.h 1856 2010-03-22 15:55:26Z maufder $

#ifndef _MDBOUNDARYMODEL_H_
#define _MDBOUNDARYMODEL_H_

#include "PhysicalModel.h"
#include "SimulationOptions.h"
#include "MaterialBoundary.h"

class Elem;
class Point;


//! The base class for MD boundary conditions
class MDBoundaryModel : public PhysicalModel
{

  public:

    //! Destructor
    ~MDBoundaryModel(void) {};

    //! Creator function
       static MDBoundaryModel* create(const MaterialBoundary* boundary,
           const ModelOptions& options);


    //! Calculate for a point on the given side
    virtual void calculate(const Elem* elem, unsigned int side,
        const Point& point) = 0;

    //! Get the depth of the region where the temperature is applied
    double get_depth(void) const;

    //! Get the temperature
    double get_temperature(void) const;

    //! Get the temperature damp
    double get_temp_damp(void) const;

    //! Get the drag factor
    double get_drag_factor(void) const;


  protected:

    //! Constructor
    MDBoundaryModel(const ModelOptions& options);

    //! Set the temperature
    void set_temperature(double temperature);

    //! Set the depth of the region where the temperature is applied
    void set_depth(double depth);

    //! Set the temperature damp
    void set_temp_damp(double Tdamp);

    //! Set the drag factor
    void set_drag_factor(double drag_factor);

  private:

    //! The temperature
    double _temperature;

    //! The depth of the region where the temperature is applied
    double _depth;

    //! The temperature damp
    double _Tdamp;

    //! The drag factor
    double _Dfactor;


};


inline
MDBoundaryModel::MDBoundaryModel(const ModelOptions& options) :
  PhysicalModel(options),
  _temperature(SimulationOptions::temperature),
  _depth(0.5),
  _Tdamp(10.0),
  _Dfactor(0.2)
{
}


inline
double
MDBoundaryModel::get_depth() const
{
  return _depth;
}

inline
double
MDBoundaryModel::get_temperature() const
{
  return _temperature;
}


inline
double
MDBoundaryModel::get_temp_damp(void) const
{
     return _Tdamp;
}

inline
double
MDBoundaryModel::get_drag_factor(void) const
{
     return _Dfactor;

}

inline
void
MDBoundaryModel::set_temperature(double temperature)
{
  _temperature = temperature;
}


inline
void
MDBoundaryModel::set_drag_factor(double Dfactor)
{
  _Dfactor = Dfactor;
}


inline
void
MDBoundaryModel::set_temp_damp(double Tdamp)
{
  _Tdamp = Tdamp;
}


inline
void
MDBoundaryModel::set_depth(double depth)
{
  _depth = depth;
}



#endif // _MDBOUNDARYMODEL_H_
