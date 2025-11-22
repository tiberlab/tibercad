// $Id$

#ifndef _EXTERNALPROFILE_H_
#define _EXTERNALPROFILE_H_

#include "tibercad/module/SolutionProvider.h"

#include "tibercad/module/TiberModelObject.h"
#include "tibercad/base/libMeshDefs.h"

class Elem;

/*!
 * \brief Base class for reading external data profiles
 */
class ExternalProfile : protected TiberModelObject
{
  public:

    //! Destructor
    virtual ~ExternalProfile(void);

    //! The creation method
    static ExternalProfile* create(const ModelOptions& options);

    //! Get the data at an element
    virtual double get_data(const Elem* elem) const;

    //! Get the data at a coordinate
    virtual double get_data(const Elem* elem, const Point& p) const;

    //! Get extremal values
    virtual std::pair<double, double> get_min_max(void) const;


  protected:

    //! Constructor
    ExternalProfile(const ModelOptions& options);


  private:

    //! Setup, called only if now derived class is instantiated
    void setup(void);

    //! The data provider
    SolutionProvider _data_source;

};



#endif //_EXTERNALPROFILE_H_
