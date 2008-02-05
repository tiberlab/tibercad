#ifndef _POISSON_H_
#define _POISSON_H_


//------------------------------------------------------------------------------

#include "SimulationInterface.h"
#include "elem.h"
#include "PoissonModel.h"
#include "linear_implicit_system.h"
#include "point.h"

class Device;
class Mesh;
//!  Class to solve Poisson problem
class Poisson : public SimulationInterface
/*!
  
 \f$ \nabla_i \cdot(<epsilon \nabla \varphi + \mathbf{P}) =-frac{\rho}{\epsilon}\f$

 * The get_solution() methods can provide the following variables:
 * Electrostatic potential (V)
*/

{
 public:
  
  //! internal variables
  enum Variables
    {
      UNKNOW,
      POTENTIAL
    };
  
  /*! \copydoc SimulationInterface::convert_variable_name_to_id() */
  virtual ID convert_variable_name_to_id(const std::string& variable_name) const;

  //!options that we need for this simulation
  struct options
  {

    double work_units; //!< SI units, has to be consistent with the database parameters

    double length_scale; //!< mesh_units/work_units

  };

  //Init the poisson Model
  void init_poisson_model(const Elem* elem);

  //!Get a avarage temperature for a given element  
  double get_potential_element(const Elem* elem) const;

  //!Get a temperature for all nodes of a given element
  std::vector<double> get_potential_node(const Elem* elem);

  //!Constructor
   Poisson(void);
  
  //!Destructor
  virtual ~Poisson(void);
  
  virtual PhysicalModel* create_physical_model(const ModelOptions &options,
      const Material* mat) const throw (ModelErrorException);
  
  
  virtual BoundaryProperties* create_boundary_model(const ModelOptions &options) const 
    throw (ModelErrorException);

  static void assemble_poisson_matrix(EquationSystems& es,
				     const std::string& system_name);

  //!Create an Poisson object 
  static  Poisson*  create(void);

  //!Get a solution given the element and the point
  void get_solution(const Elem* elem, const std::vector<Point>& p,
		    std::vector<double>& solution);
  
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


  
 private:
  
  PoissonModel* poisson_model; 
    

   //! Quadrature point along the face of the element 
  const std::vector<Point> qface_point;


  EquationSystems* 	equation_systems;

  std::string system_name;
  
  LinearImplicitSystem* my_system;  
 
  //! Order the solution in correct mode
  virtual void 	build_nodal_results(const std::set< std::string > &variables, 
				     std::vector< double > &results, 
				     std::vector< std::string > &legend);

  
  static Poisson* static_this;

  options opt;

  //!non-static method that actually does matrix assembling 
  void do_assemble(EquationSystems& es, const std::string& system_name);

  static Device* _device;

  //!Dimension of meshmap
  short dim;  

  //!Pointer to mesh
  Mesh* mesh;

 protected:
  
 
  //! \copydoc  SimulationInterface::do_init() 
  virtual void 	do_init (void);
  
  //!Do the solve
  virtual void do_solve (void);
 
  //!Parse the options
  virtual void 	parse_options(void);
 

};
 
#endif
