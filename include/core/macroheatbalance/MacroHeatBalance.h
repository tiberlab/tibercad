// $Id$

#ifndef _MACROHEATBALANCE_H_
#define _MACROHEATBALANCE_H_



//------------------------------------------------------------------------------


#include "SimulationInterface.h"

#include "HeatModel.h"


class TiberLinearSystem;
class Device;
class Mesh;


//!  Class to solve heat transport problem
class MacroHeatBalance : public SimulationInterface
/*!
  
 \f$ -\nabla_i \cdot(k_{i,j} \nabla_j T)+ \nabla \cdot (T P_n J_n + T P_p J_p)=-\nabla cdot (J_n \phi_n + J_p \phi_p)\f$
 * The get_solution() methods can provide the following variables:
 * Temperature (K)
*/

{

 public:


   enum Variables
   {
     UNKNOWN = 0,
     TEMPERATURE,
     JQX,
     JQY,
     JQZ
   };
       




  //!options that we need for this simulation
  struct options
  {
  
    std::string  kappa_solve; //!< Model for lattice thermal conductivity
      
    double max_error; //!< Max tollerance for self-consistent loop  

    double work_units; //!< SI units, has to be consistent with the database parameters

    /**
     * The order of gauss integration
     */
    libMeshEnums::Order integration_order;

  };

    
  /*!
   * \copydoc SimulationInterface::get_solution_secure(const Elem*,
   * const std::vector<Point>&, const std::vector<ID>&,
   * std::vector<std::vector<double> >&)
   */
  virtual void get_solution_secure(const Elem* elem, const std::vector<Point>& p,
				   const std::set<ID>& ids, std::vector<std::map<ID, double> >& values);
  
  
  /*!
   * \copydoc SimulationInterface::get_solution_secure(const Elem*,
   * const std::vector<ID>&, std::vector<std::vector<double> >&)
   */  
  virtual void get_solution_secure(const Elem* elem,
				   const std::set<ID>& ids, std::vector<std::map<ID, double> >& values);
  

  //!Constructor
  MacroHeatBalance();
  
  //!Destructor
  virtual ~MacroHeatBalance();
  
  virtual PhysicalModel* create_physical_model(const ModelOptions &options,
      const Material* mat) const throw (ModelErrorException);
  
  
  virtual BoundaryProperties* create_boundary_model(const ModelOptions &options) const 
    throw (ModelErrorException);

  static void assemble_heat_matrix(EquationSystems& es,
				     const std::string& system_name);

  //!Create an MacroHeatBalance object 
  static MacroHeatBalance*  create(void);
  
  /*! \copydoc SimulationInterface::build_integrated_quantities() */
  virtual void build_integrated_quantities(
					   const std::set<std::string>& names,
					   std::vector<double>& values);
  
  
  /*! \copydoc SimulationInterface::build_integrated_quantities_description()
   */
  virtual void build_integrated_quantities_description(
						       const std::set<std::string>& names,
						       std::vector<std::string>& legend,
						       std::vector<std::string>& description);

   /*! \copydoc SimulationInterface::do_print_info() */
    virtual void do_print_info(void);

 private:
  
  std::string heat_legend;

  std::set<ID> JQ_var;

   //! Quadrature point along the face of the element 
  const std::vector<Point> qface_point;

  //! A pointer to heat simulation
  MacroHeatBalance* _heat_simul;

  EquationSystems * 	equation_systems;
  
  HeatModel* heat_model; 

  //void init_heat_model(const Elem* elem);
 

  std::string system_name;
  
  TiberLinearSystem* my_system;  
 
  //! Order the solution in correct mode
  virtual void 	build_nodal_results(const std::set< std::string > &variables, 
				     std::vector< double > &results, 
				     std::vector< std::string > &legend);

   
  //! Order the solution in correct mode
  virtual void build_elemental_results(const std::set<std::string>& variables,
				       std::vector<double>& results, 
				       std::vector<std::string>& legend);

  //! Calculate Power Dissipated 
  /*!
   * Integrates numerically over the boundary elements.
   * The power dissipated is then:
   *
   * \f[P = \int_{\Gamma} -\kappa \nabla T \cdot \mathbf{N} \mathrm{d}\Gamma \f]
   */
  void calculate_power_surfint(void);
  
  static MacroHeatBalance* static_this;

  options myopts;

  //!non-static method that actually does matrix assembling 
  void do_assemble(EquationSystems& es, const std::string& system_name);

  static Device* _device;

  //!Dimension of meshmap
  short dim;  

  //!Power Dissipated
  double _power;

  //!Pointer to mesh
  Mesh* mesh;

 protected:

    /*! \copydoc SimulationInterface::convert_variable_name_to_id() */
  virtual ID convert_variable_name_to_id(const std::string& variable_name) const;
 
  //! \copydoc  SimulationInterface::do_init() 
  virtual void 	do_init (void);
  
  //!Do the solve
  virtual void do_solve (void);
 
  //!Parse the options
  virtual void 	parse_options(void);
 

};
 
#endif
