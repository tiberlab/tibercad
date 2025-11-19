#ifndef _DEGRADATIONMODEL_H_
#define _DEGRADATIONMODEL_H_

#include "tibercad/physics/PhysicalModel.h"

/*!
 * \brief Base class for degradation models
 *
 * The model is intended for use with solar cell
 * equivalent circuits. It has to implement a method
 * that returns the model parameters changed according
 * to some degradation model.
 */
class DegradationModel : public PhysicalModel
{

  public:

    //! A base class for parameter containers
    /*!
     * Probably double values is what is needed, but
     * derived classes might need other stuff.
     */
    class Parameters {
      public:
        std::vector<double> double_params;
    };

    virtual ~DegradationModel(void) = default;

    //! Implementation for the 1-diode equivalent circuit
    void degrade_parameters(const libMesh::Elem* elem,
                            const libMesh::Point& p,
                            Parameters& params) const
    {
      do_degrade_params(elem, p, params);
    };


  protected:

    //! Constructor
    DegradationModel(const ModelOptions& options)
      : PhysicalModel(options) {};


    virtual void do_degrade_params(const libMesh::Elem* elem,
                                   const libMesh::Point& p,
                                   Parameters& params) const = 0;

  private:

};


#endif // _DEGRADATIONMODEL_H_
