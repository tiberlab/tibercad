
#ifndef _CONSTANTDOS_H_
#define _CONSTANTDOS_H_


#include "tibercad/model_base/DensityOfStates.h"

class TBDLEXPORT ConstantDOS : public DensityOfStates
{

  public:

    //! Destructor
    virtual ~ConstantDOS(void) {};

    //! Creator function
    static ConstantDOS* create(const ModelOptions& options);

  protected:

    //! Constructor
    ConstantDOS(const ModelOptions& options);

    virtual void read_database(void);

    //! Get occupied states and the derivative with respect to phi
    virtual void
    calculate_density_and_derivative(std::vector<double>& result, double E, double Epot,
                                      double kT, double kTlattice, const Elem* elem, const Point& p) const;

    virtual void do_init(void);
	
  private:
  	
    //
    double _Ewidth;

    //Total density parameter
    double _N0;
};

//
// inline methods
//

inline
ConstantDOS*
ConstantDOS::create(const ModelOptions& options)
{
  return new ConstantDOS(options);
}


#endif // _CONSTANTDOS_H_
