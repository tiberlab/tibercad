// $Id$

#ifndef _COMPOSITEROFILE_H_
#define _COMPOSITEROFILE_H_


#include "tibercad/profiles/ExternalProfile.h"

#include "point.h"
#include "elem.h"

#include <vector>

//! Class to create profiles in multiple directions
class CompositeProfile : public ExternalProfile
{
  public:

    //! Constructor
    CompositeProfile(const ModelOptions& options);

    //! Destructor
    virtual ~CompositeProfile(void);

    virtual double get_data(const libMesh::Elem* elem, const libMesh::Point& p) const override;

    virtual std::pair<double, double> get_min_max(void) const override;

  private:

    //! The peak value, used as multiplier
    double _peak;

    //! An offset, will be added to the multiplied value
    double _offset;

    //! The list of the profiles
    std::vector<ExternalProfile*> _profiles;
};


#endif //_COMPOSITEPROFILE_H_
