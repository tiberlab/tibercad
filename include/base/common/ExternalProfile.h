// $Id$

#include "TiberModelObject.h"

class Point;

/*!
 * \brief Base class for reading external data profiles
 */
class ExternalProfile : protected TiberModelObject
{
  public:

    //! Destructor
    virtual ~ExternalProfile(void);

    //! The creation method
    ExternalProfile* create(const ModelOptions& options);

    //! Get the data at a coordinate
    virtual double get_data(const Point& p) const = 0;

  private:

    //! Constructor
    ExternalProfile(const ModelOptions& options);

};
