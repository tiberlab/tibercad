// $Id: ExponentialDOS.h 2117 2010-10-20 20:56:13Z maufder $

#ifndef _EXPONENTIALDOS_H_
#define _EXPONENTIALDOS_H_


#include "DensityOfStates.h"


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



    //! Get occupied states
    /*!
     * \return the density of occupied states in cm^-3
     * for the given energy \c E, where \c E is something
     * like \$E_0 - E_f\$.
     */
    virtual double get_occupied_density(double E, double kT) const;


    //! Get the derivative with respect to the argument
    virtual double get_occupied_density_derivative(double E, double kT) const;


  protected:

    //! Constructor
    ExponentialDOS(const ModelOptions& options);

    virtual void do_init(void);


  private:

    //! The tail parameter
    double _alpha;

    //! Calculates the value under the integral
    double get_value(double e, double E, double kT) const;

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
