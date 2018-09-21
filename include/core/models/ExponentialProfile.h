// $Id$

#ifndef _EXPONENTIALPROFILE_H_
#define _EXPONENTIALPROFILE_H_


#include "ExternalProfile.h"

#include "point.h"
#include "elem.h"

#include <vector>

//! Class to create a linear profile
class ExponentialProfile : public ExternalProfile
{
  public:

    //! Constructor
    ExponentialProfile(const ModelOptions& options);

    //! Destructor
    virtual ~ExponentialProfile(void);

    virtual double get_data(const Elem* elem, const Point& p) const override;

    virtual std::pair<double, double> get_min_max(void) const override;

  private:

    double _max;

    double _lambda;

    enum Type {
      onesided,
      symmetric,
      continued,
    };

    Type _type;

    Point _direction;
    Point _origin;

};


#endif //_EXPONENTIALPROFILE_H_
