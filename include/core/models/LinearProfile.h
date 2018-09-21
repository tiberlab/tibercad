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

    virtual double get_data(const libMesh::Elem* elem, const libMesh::Point& p) const override;

    virtual std::pair<double, double> get_min_max(void) const override;

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

    libMesh::Point _direction;
    libMesh::Point _origin;

};


#endif //_LINEARPROFILE_H_
