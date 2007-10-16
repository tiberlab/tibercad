#ifndef _DFTB_H_
#define _DFTB_H_

using namespace std;

//---------------------------------------------------------------------

#include "TightBinding.h"
#include "AtomisticStructure.h"
#include "SimulationEnvironment.h"


class Dftb : public TightBinding{

public:

struct dftb_options
  {
    vector<string> sk_files;
  };


 //! Constructor
  Dftb(void);

  //! Destructor
  ~Dftb(void);

  //! Create TightBinding object
  static Dftb* create();


  private:

  //! Get options suited for DFTB+ tight binding builder and solver
  void get_dftbp_options();

  //! Structure containing options for DFTB+ tight binding builder and solver
  dftb_options _dftb_options;

  //! Pointer to atomistic structure for the simulation
  AtomisticStructure* _atomistic_structure;

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
