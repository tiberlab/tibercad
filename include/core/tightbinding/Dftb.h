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

struct dftb_options
  {
    std::vector<std::string> sk_files;
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

  //! A function for building a DFTB compatible char of SK names,
  //! based on species
  char* build_sk_names(void);

 //! Structure containing options for DFTB+ tight binding builder and solver
  dftb_options _dftb_options;


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
