// $Id: BulkDOS.h 2117 2010-10-20 20:56:13Z maufder $

#ifndef _BULKDOS_H_
#define _BULKDOS_H_


#include "tibercad/physics/misc/DensityOfStates.h"


/*!
 * \brief Base class for density of states
 */
class TBDLEXPORT BulkDOS : public DensityOfStates
{

  public:

    //! Destructor
    virtual ~BulkDOS(void) {};


    //! Creator function
    static BulkDOS* create(const ModelOptions& options);



  protected:

    //! Constructor
    BulkDOS(const ModelOptions& options);

    //! Read band edge, mass, degeneracy from database
    /*!
     * This will only read from the database if the \c particle option
     * is provided as \c electron or \c hole
     */
    virtual void read_database(void);

    virtual void do_init(void);

    virtual void do_print_info(void);


    //! Get occupied states and derivative w.r.t. Ef
    virtual void
    calculate_density_and_derivative(std::vector<double>& result, double Ef, double Epot,
        double kT, double kTlattice, const Elem* elem, const Point& p) const;


  private:

    //! A structure for band gap parameters
    struct GapParameters
    {
      double Ev;
      double Eg0;
      double varshni_a;
      double varshni_b;
      void zero(void) { Ev = Eg0 = varshni_a = varshni_b = 0.0; }
      double gap(double T) const { return Eg0 - varshni_a * T * T / (T + varshni_b); }
    };


    //! We allow for several subbands
    std::vector<double> _ref_energies;

    //! The effective masses
    std::vector<double> _dos_mass;

    //! The degeneracies
    std::vector<int> _degeneracy;

    //! The band gap parameters
    std::vector<GapParameters> _gap_params;


    //! The DOS factor
    /*!
     * In 3D the effective DOS is
     * \f[ N_{3D} = \f]
     */
    double _dos_factor;


};

//
// inline methods
//

inline
BulkDOS*
BulkDOS::create(const ModelOptions& options)
{
  return new BulkDOS(options);
}


#endif // _BULKDOS_H_
