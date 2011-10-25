// $Id: DensityOfStates.h 2117 2010-10-20 20:56:13Z maufder $

#ifndef _DENSITYOFSTATES_H_
#define _DENSITYOFSTATES_H_


#include "PhysicalModelInterface.h"


/*!
 * \brief Base class for density of states
 */
class TBDLEXPORT DensityOfStates : public PhysicalModelInterface
{

  public:

    //! Destructor
    virtual ~DensityOfStates(void) {};


    //! Creator function
    static DensityOfStates* create(const ModelOptions& options);



    //! Get occupied states
    /*!
     * \return the density of occupied states in cm^-3
     * for the given energy \c E, where \c E is something
     * like \$E_0 - E_f\$.
     */
    virtual double get_occupied_density(double E, double kT) const = 0;


    //! Get the derivative with respect to the argument
    virtual double get_occupied_density_derivative(double E, double kT) const = 0;


  protected:

    //! Constructor
    DensityOfStates(const ModelOptions& options);


  private:


};

//
// inline methods
//



#endif // _DENSITYOFSTATES_H_
