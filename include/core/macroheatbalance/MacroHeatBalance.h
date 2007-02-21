#ifndef _MACROHEATBALANCE_H_
#define _MACROHEATBALANCE_H_



//------------------------------------------------------------------------------


#include "SimulationInterface.h"


class DriftDiffusion;

//!  Class to solve heat transport problem
class MacroHeatBalance : public SimulationInterface
/*!
   Class to solve heat transport problem
*/

{
 public:
  //!Constructor
  MacroHeatBalance();
  
  //!Destructor
  ~MacroHeatBalance();
  
  PhysicalModel* create_physical_model(const ModelOptions &options) const 
    throw (ModelErrorException);
  
  
  BoundaryProperties* create_boundary_model(const ModelOptions &options) const 
    throw (ModelErrorException);

  
  
  static MacroHeatBalance*  create(void);
  
  
 
 private:
  
  virtual void 	build_nodal_results(const std::set< std::string > &variables, 
				     std::vector< double > &results, 
				     std::vector< std::string > &legend);

  //! Pointer to a DriftDiffusion simulation
  DriftDiffusion* _dd_simul;

 protected:
  
  virtual void 	do_init (void);
  
  virtual void do_solve (void);
 
  virtual void 	parse_options(void);
 

};
 
#endif
