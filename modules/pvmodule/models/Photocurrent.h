#ifndef _PHOTOCURRENT_H_
#define _PHOTOCURRENT_H_

#include "tibercad/physics/PhysicalModel.h"

/*!
 * \brief Base class for photocurrent models
 *
 * The model has to implement a method that
 * returns the photocurrent density in function
 * of a given coordinate.
 */
class Photocurrent : public PhysicalModel
{

  public:

    virtual ~Photocurrent(void) {};

    double get_photocurrent(const libMesh::Elem* elem,
                            const libMesh::Point& p) const
    {
      return do_get_photocurrent(elem, p);
    };


  protected:

    //! Constructor
    Photocurrent(const ModelOptions& options)
      : PhysicalModel(options) {};


    virtual double do_get_photocurrent(const libMesh::Elem* elem,
                                       const libMesh::Point& p) const = 0;

  private:

};


#endif // _PHOTOCURRENT_H_
