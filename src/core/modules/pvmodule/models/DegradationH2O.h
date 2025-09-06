#ifndef _DEGRADATIONH2O_H_
#define _DEGRADATIONH2O_H_

#include "DegradationModel.h"
#include "SolutionProvider.h"

/*!
 * \brief An example for degradation due to relative humidity
 *
 * This class implements an equiv. circuit parameter dependency 
 * on relative humidity, calculated from a water ingress model.
 * The formulas have been obtained by fitting to data in
 * Bhatt et al., Organic Electronics 39 (2016) 258e266.
 */
class DegradationH2O : public DegradationModel
{

  public:

    virtual ~DegradationH2O(void) = default;

    static DegradationH2O* create(const ModelOptions& options);


  protected:

    virtual void do_init(void) override;

    virtual void do_degrade_params(const libMesh::Elem* elem,
                                   const libMesh::Point& p,
                                   DegradationModel::Parameters& params) const final;


  private:

    //! Private constructor
    /*!
     * This is a specific model. Other models should directly
     * derive from the base class to not create a mess
     */
    DegradationH2O(const ModelOptions& options);

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


#endif // _DEGRADATIONH2O_H_
