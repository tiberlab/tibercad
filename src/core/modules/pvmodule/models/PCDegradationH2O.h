#ifndef _PCDEGRADATIONH2O_H_
#define _PCDEGRADATIONH2O_H_

#include "Photocurrent.h"
#include "SolutionProvider.h"

/*!
 * \brief An example for photocurrent dependency on H20 concentration
 *
 * This class implements a photocurrent dependency based on the
 * assumption, that presence of water in a perovskite cell leads
 * to material degradation and reduction of optically active
 * material.
 */
class PCDegradationH2O : public Photocurrent
{

  public:

    ~PCDegradationH2O(void) {};

    static PCDegradationH2O* create(const ModelOptions& options);


  protected:

    virtual void do_init(void) final;

    virtual double do_get_photocurrent(const libMesh::Elem* elem,
                                       const libMesh::Point& p) const final;


  private:

    PCDegradationH2O(const ModelOptions& options);

    //! The current degradation factor
    /*!
     * The formula is based on a fit of the data in [ref].
     * It uses a generalized logistic function.
     */
    double degradation_factor(double humidity) const;

    //! The initial (undegraded) photocurrent
    double _initial_current = 0.02;

    //! The reference humidity in the degradation fit
    double _RH_ref = 72.1;

    //! The exponent in the degradation fit
    double _exponent = 8.28;

    //! From where to get relative humidity
    SolutionProvider _humidity_model;

};


#endif // _PCDEGRADATIONH2O_H_
