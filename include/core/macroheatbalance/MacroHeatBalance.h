#ifndef _MACROHEATBALANCE_H_
#define _MACROHEATBALANCE_H_



//------------------------------------------------------------------------------


#include "SimulationInterface.h"

#include "linear_implicit_system.h"
class DriftDiffusion;
class Device;
class Mesh;
//!  Class to solve heat transport problem
class MacroHeatBalance : public SimulationInterface
/*!
   
*/

{
 public:

  //!options that we need for this simulation
  struct options
  {
    double  lin_tol; //!< linear tolerance 
    std::string  current_simulation; //!< name of drift-diffusion simmulation
    bool Thompson_effect; 
    
    double work_units; //!< SI units, has to be consistent with the database parameters

    double length_scale; //!< mesh_units/work_units

  };

  //!Constructor
  MacroHeatBalance();
  
  //!Destructor
  virtual ~MacroHeatBalance();
  
  virtual PhysicalModel* create_physical_model(const ModelOptions &options) const 
    throw (ModelErrorException);
  
  
  virtual BoundaryProperties* create_boundary_model(const ModelOptions &options) const 
    throw (ModelErrorException);

  static void assemble_heat_matrix(EquationSystems& es,
				     const std::string& system_name);

  //!Create an MacroHeatBalance object 
  static MacroHeatBalance*  create(void);
  
  
 private:
  
  EquationSystems * 	equation_systems;
  

  std::string system_name;

  LinearImplicitSystem* my_system;  
 
  
  virtual void 	build_nodal_results(const std::set< std::string > &variables, 
				     std::vector< double > &results, 
				     std::vector< std::string > &legend);

  //! Pointer to a DriftDiffusion simulation
  DriftDiffusion* _dd_simul;

  
  static MacroHeatBalance* static_this;

  options opt;

  //!non-static method that actually does matrix assembling 
  void do_assemble(EquationSystems& es, const std::string& system_name);

  
  static Device* _device;

  //!Dimension of mesh
  short dim;  

  //!Pointer to mesh
  Mesh* mesh;

 protected:
  
 
  //! \copydoc  SimulationInterface::do_init() 
  /*!
    In this class it does ...
   */
  virtual void 	do_init (void);
  
  //!Do the solve
  virtual void do_solve (void);
 
  //!Parse the options
  virtual void 	parse_options(void);
 

};
 
#endif
