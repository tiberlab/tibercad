// $Id: ExponentialDOS.h 2117 2010-10-20 20:56:13Z maufder $

#ifndef _EXPONENTIALDOS_H_
#define _EXPONENTIALDOS_H_


#include "tibercad/model_base/DensityOfStates.h"


/*!
 * \brief Base class for density of states
 */
class TBDLEXPORT ExponentialDOS : public DensityOfStates
{

  public:

    //! Destructor
    virtual ~ExponentialDOS(void) {};


    //! Creator function
    static ExponentialDOS* create(const ModelOptions& options);



  protected:

    //! Constructor
    ExponentialDOS(const ModelOptions& options);

    virtual void do_init(void);


    //! Get occupied states and derivative
    virtual void
    calculate_density_and_derivative(std::vector<double>& result, double E, double Epot,
        double kT, double kTlattice, const Elem* elem, const Point& p) const;



  private:

    //! The tail parameter
    double _alpha;

    //! Calculates the value under the integral
    double _get_value(double e, double E, double kT) const;


};

//
// inline methods
//

inline
ExponentialDOS*
ExponentialDOS::create(const ModelOptions& options)
{
  return new ExponentialDOS(options);
}


#endif // _EXPONENTIALDOS_H_
