// $Id$

#ifndef _EXTPROFILE1D_H_
#define _EXTPROFILE1D_H_


#include "ExternalProfile.h"

#include "point.h"

#include <vector>

//! Class to read an external line profile
class ExtProfile1D : public ExternalProfile
{
  public:

    //! Constructor
    ExtProfile1D(const ModelOptions& options);

    //! Destructor
    virtual ~ExtProfile1D(void);

    virtual double get_data(const Point& p) const;

    virtual std::pair<double, double> get_min_max(void) const;

  private:

    void _read_source(void);

    std::vector<double> _x_coord;
    std::vector<double> _values;
    std::vector<size_t> _addressing;

    double _min;
    double _max;

    Point _direction;
    Point _origin;
    double _scale;
    double _data_scale;

};


#endif //_EXTPROFILE1D_H_
