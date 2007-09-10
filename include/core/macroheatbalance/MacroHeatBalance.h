#ifndef _MACROHEATBALANCE_H_
#define _MACROHEATBALANCE_H_



//------------------------------------------------------------------------------


#include "SimulationInterface.h"

#include "HeatModel.h"
#include "linear_implicit_system.h"
class DriftDiffusion;
class Device;
class Mesh;
//!  Class to solve heat transport problem
class MacroHeatBalance : public SimulationInterface
/*!
  
 \f$ -\nabla_i \cdot(k_{i,j} \nabla_j T)+ \nabla \cdot (T P_n J_n + T P_p J_p)=-\nabla cdot (J_n \phi_n + J_p \phi_p)\f$
*/

{
 public:


 enum Variables
     {
            UNKNOWN = INVALID_ID,
            TEMPERATURE
     };
       




  //!options that we need for this simulation
  struct options
  {
    double  lin_tol; //!< linear tolerance 
    std::string  current_simulation; //!< name of drift-diffusion simmulation
    std::string  thomson_peltier_effect;
  
    std::string  kappa_solve; //!< Model for lattice thermal conductivity
      
    double max_error; //!< Max tollerance for self-consistent loop  

    double work_units; //!< SI units, has to be consistent with the database parameters

    double length_scale; //!< mesh_units/work_units

  };

  
  
  virtual void get_solution_secure(const Elem* elem, const std::vector<Point>& p,
				   const std::set<ID>& ids, std::vector<std::map<ID, double> >& values);
   

   
virtual void get_solution_secure(const Elem* elem,
				 const std::set<ID>& ids, std::vector<std::map<ID, double> >& values);


  //!Get a avarage temperature for a given element  
  double get_temperature_element(const Elem* elem) const;

  //!Get a temperature for all nodes of a given element
  std::vector<double> get_temperature_node(const Elem* elem);

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
  
 


   //! Quadrature point along the face of the element 
  const std::vector<Point> qface_point;

  //! A pointer to heat simulation
  MacroHeatBalance* _heat_simul;

  EquationSystems * 	equation_systems;
  
  HeatModel* heat_model; 

  void init_heat_model(const Elem* elem);
 

  std::string system_name;
  
  LinearImplicitSystem* my_system;  
 
  //! Order the solution in correct mode
  virtual void 	build_nodal_results(const std::set< std::string > &variables, 
				     std::vector< double > &results, 
				     std::vector< std::string > &legend);

  
  static MacroHeatBalance* static_this;

  options opt;

  //!non-static method that actually does matrix assembling 
  void do_assemble(EquationSystems& es, const std::string& system_name);

  static Device* _device;

  //!Dimension of meshmap
  short dim;  

  //!Pointer to mesh
  Mesh* mesh;

 protected:

  //! The correct keyword is "temperature"
  virtual ID convert_variable_name_to_id(const std::string& variable_name);
 
  //! \copydoc  SimulationInterface::do_init() 
  virtual void 	do_init (void);
  
  //!Do the solve
  virtual void do_solve (void);
 
  //!Parse the options
  virtual void 	parse_options(void);
 

};
 
#endif
