// $Id$

#ifndef _LINEARPROFILE_H_
#define _LINEARPROFILE_H_


#include "ExternalProfile.h"

#include "point.h"
#include "elem.h"

#include <vector>

//! Class to create a linear profile
class LinearProfile : public ExternalProfile
{
  public:

    //! Constructor
    LinearProfile(const ModelOptions& options);

    //! Destructor
    virtual ~LinearProfile(void);

    virtual double get_data(const Elem* elem, const Point& p) const;

    virtual std::pair<double, double> get_min_max(void) const;

  private:

    double _min;
    double _max;

    double _distance;

    enum Type {
      onesided,
      symmetric,
      continued,
    };

    Type _type;

    Point _direction;
    Point _origin;

};


#endif //_LINEARPROFILE_H_
