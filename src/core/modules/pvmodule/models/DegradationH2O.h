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
 * 
 * The fit functions are power laws, in particular:
 * 
 * \f{eqnarray*}
 *  \frac{R_s}{R_{s,0}} & = & 1 + \left(\frac{RH}{RH_0})^\gamma \\
 *  \frac{I_ph}{I_{ph,0}} & = & \frac{1}{1 + \left(\frac{RH}{RH_0})^\gamma} \\
 * \f}
 * 
 * Here \f$ RH \f$ is the relative humidity, which has to be provided by
 * another module.
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

    //! The reference humidity in the photocurrent degradation fit
    double _RH_ref_ph = 1e9;

    //! The exponent in the photocurrent degradation fit
    double _exponent_ph = 1;

    //! The reference humidity in the series resistance degradation fit
    double _RH_ref_rs = 1e9;

    //! The exponent in the series resistance degradation fit
    double _exponent_rs = 1;

    //! From where to get relative humidity
    SolutionProvider _humidity_model;

};


#endif // _DEGRADATIONH2O_H_
