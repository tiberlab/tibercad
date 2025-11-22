// $Id$

#ifndef _DELTADOS_H_
#define _DELTADOS_H_


#include "tibercad/physics/misc/DensityOfStates.h"


/*!
 * \brief Base class for density of states
 */
class TBDLEXPORT DeltaDOS : public DensityOfStates
{

  public:

    //! Destructor
    virtual ~DeltaDOS(void) {};


    //! Creator function
    static DeltaDOS* create(const ModelOptions& options);




  protected:

    //! Constructor
    DeltaDOS(const ModelOptions& options);

    virtual void do_init(void);


    //! Get occupied states and derivative w.r.t. Ef
    virtual void
    calculate_density_and_derivative(std::vector<double>& result, double Ef, double Epot,
        double kT, double kTlattice, const Elem* elem, const Point& p) const;


  private:

  //! The effective DOS
  double _N0;
};

//
// inline methods
//

inline
DeltaDOS*
DeltaDOS::create(const ModelOptions& options)
{
  return new DeltaDOS(options);
}


#endif // _DELTADOS_H_
