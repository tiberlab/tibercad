#ifndef _DFTB_H_
#define _DFTB_H_

//---------------------------------------------------------------------

#include "TightBinding.h"
#include "AtomisticStructure.h"
#include "SimulationEnvironment.h"

//! A class for Tight Binding simulations with DFTB+ code
/*!
 *This class provides methods for invoking DFTB+ (Density Functional
 *Tight Binding) library, in order to make calculations at equilibrium
 */
class Dftb : public TightBinding{

public:

struct DftbOptions
  {
    char* skNames;
    char* speciesNames;
    int nAtom;
    int nType;
    double* coords;
    int* species;
    int iPeriodic;

  };

 //! Constructor
  Dftb(void);

  //! Destructor
  ~Dftb(void);

  //! Create object
  static Dftb* create();


  private:


  //! Get options suited for DFTB+ tight binding builder and solver
  void get_dftbp_options();

  //! A function for building a DFTB compatible char of SK and species names
  void build_names(void);

  //! Function for building options from atomistic structure
  void build_structure_options(void);

 //! Structure containing options for DFTB+ tight binding builder and solver
  DftbOptions _dftb_options;



protected:

  virtual void do_init(void);

  virtual void do_solve (void);

  virtual void  parse_options(void);

};



inline 
Dftb* Dftb::create()
{
  return new Dftb;
}





#endif
