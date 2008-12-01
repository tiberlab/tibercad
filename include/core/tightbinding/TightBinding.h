#ifndef _TIGHTBINDING_H_
#define _TIGHTBINDING_H_

//-----------------------------------------------------------------------------------------

#include "SimulationInterface.h"
#include "AtomisticStructure.h"


//forward declaration
class Device;
class Mesh;
class DftbpWrapper;

//!Main class for Atomistic Tight Binding simulation at equilibrium
//! DFTB code is used for simulations
class TightBinding : public SimulationInterface{


public:

	enum Shell
	{
		NONE = 0,
		S = 1,
		P = 2,
		D = 3
	};

  //! Constructor
  TightBinding();

  //! Destructor
  ~TightBinding();

  //! Create TightBinding object
  static TightBinding* create();

  virtual PhysicalModel* create_physical_model(const ModelOptions &options,
      const Material* mat) const throw (ModelErrorException);

  virtual BoundaryProperties* create_boundary_model(const ModelOptions &options) const
    throw (ModelErrorException);



private:

protected:

  virtual void  do_init (void);

  virtual void do_solve (void);

  virtual void  parse_options(void);

  virtual void obtain_hubbard_parameters(void);

  //! Pointer to atomistic structure for the simulation;
  AtomisticStructure* _atomistic_structure;

  //! Get the atomistic structure pointer from the name specified in input
  //! and fill the private member _atomistic_structure
  void get_atomistic_structure(void);

//! Map of map containing hubbard parameters for any specie and any shell
/*!
 * Usage: _u_hub[<specie>][shell] = hubbard_index
 */
std::map<std::string, std::map<Shell, double> > _u_hub;

//! Build charge density on given point
double build_rho(const double x, const double y, const double z);

//! Charge variation (Mulliken Analisys) on each atom
std::vector<double> _mulliken_netcharges;

};


inline
TightBinding* TightBinding::create()
{
  return new  TightBinding();
}


#endif
