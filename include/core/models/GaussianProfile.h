// $Id$

#ifndef _GAUSSIANPROFILE_H_
#define _GAUSSIANPROFILE_H_


#include "ExternalProfile.h"

#include "point.h"
#include "elem.h"

#include <vector>

//! Class to create a Gaussian profile
class GaussianProfile : public ExternalProfile
{
  public:

    //! Constructor
    GaussianProfile(const ModelOptions& options);

    //! Destructor
    virtual ~GaussianProfile(void);

    virtual double get_data(const Elem* elem, const Point& p) const override;

    virtual std::pair<double, double> get_min_max(void) const override;

  private:

    double _max;

    double _offset;

    double _sigma;

    enum Type {
      onesided,
      symmetric,
      continued,
    };

    Type _type;

    Point _direction;
    Point _origin;

};


#endif //_GAUSSIANPROFILE_H_
