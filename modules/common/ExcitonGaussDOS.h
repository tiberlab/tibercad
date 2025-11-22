#ifndef _EXCITONGAUSSDOS_H_
#define _EXCITONGAUSSDOS_H_


#include "tibercad/physics/misc/DensityOfStates.h"


/*!
 * \brief Base class for density of states
 */
class TBDLEXPORT ExcitonGaussDOS : public DensityOfStates
{

  public:

    //! Destructor
    virtual ~ExcitonGaussDOS(void) {};


    //! Creator function
    static ExcitonGaussDOS* create(const ModelOptions& options);



  protected:

    //! Constructor
    ExcitonGaussDOS(const ModelOptions& options);

    //! Read parameters from database
    virtual void read_database(void);

    virtual void do_init(void);

    virtual void do_reinit(void);

    virtual void do_print_info(void);


    //! Get occupied states and derivative w.r.t. Ef
    virtual void
    calculate_density_and_derivative(std::vector<double>& result, double Ef, double Epot,
        double kT, double kTlattice, const Elem* elem, const Point& p) const;



  private:

    //! The exciton spin (for convenience)
    double _J;

    //! The exciton binding energy (eV)
    double _R;

    //! The gaussian DOS variance
    double _sigma;

    //! The exciton energy without considering binding energy (= gap)
    double _energy;

    //! Trapezoidal integration
    double _trapez(std::vector<double>& x, std::vector<double>& y) const;

    //! Bose - Einstein occupation function
    double _f(double x, double Ef, double kT) const;

    double _order;


};

//
// inline methods
//

inline
ExcitonGaussDOS*
ExcitonGaussDOS::create(const ModelOptions& options)
{
  return new ExcitonGaussDOS(options);
}


#endif // _BULKDOS_H_
