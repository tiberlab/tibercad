// $Id$

#ifndef _PHONONDISPERSION_H_
#define _PHONONDISPERSION_H_

#include "BoundaryProperties.h"
#include "PhysicalModel.h"
#include "mesh.h"
#include "elem.h"

//------------------------------------------------------------------------------
 

#include "SimulationInterface.h"

//#include "complex.h"

class TiberLinearSystem;
class Device;
class MeshBase;


//!  Class to solve heat transport problem
class PhononDispersion : public SimulationInterface
{

 public:


  static Device* _device;
  

  //!Pointer to mesh
  MeshBase* mesh;


    
  /*!
   * \copydoc SimulationInterface::get_solution_secure(const Elem*,
   * const std::vector<Point>&, const std::vector<ID>&,
   * std::vector<std::vector<double> >&)
   */
  virtual void get_solution_secure(const Elem* elem, const std::vector<Point>& p,
				   const std::set<ID>& ids, std::vector<std::map<ID, double> >& values){};
  
  
  /*!
   * \copydoc SimulationInterface::get_solution_secure(const Elem*,
   * const std::vector<ID>&, std::vector<std::vector<double> >&)
   */  
  virtual void get_solution_secure(const Elem* elem,
				   const std::set<ID>& ids, std::vector<std::map<ID, double> >& values){};
  

  //!Constructor
  PhononDispersion(const ModelOptions& options);
  
  //!Destructor
  virtual ~PhononDispersion();
  
   virtual PhysicalModel* create_physical_model(const ModelOptions &options,
       const Material* mat) const throw (ModelErrorException);
  
  
  virtual BoundaryProperties* create_boundary_model(const ModelOptions &options) const 
    throw (ModelErrorException);

  //!Create an MacroHeatBalance object 
  static PhononDispersion*  create(const ModelOptions& options);
  

   /*! \copydoc SimulationInterface::do_print_info() */
  virtual void do_print_info(void){};

 private:
  



  std::vector<std::vector< double > > D; //Matrix
  std::vector<double> E;                          //eignvalue
  std::vector<std::vector< double > > _eignvectors; //eignvector
  
  Tensor2 dynamical_matrix;
 
  //! Order the solution in correct mode
  virtual void 	build_nodal_results(const std::set< std::string > &variables, 
				     std::vector< double > &results, 
				    std::vector< std::string > &legend){};
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
 
  //   void diagonalize_complex(void);
    void diagonalize_double(void);
    void diagonalize_tensor(void);
    void solve_phonon_dispersion(void);
    std::map<const Elem*,std::vector<double> > PD_sol;
    std::map<const Elem*,std::vector<double> > PD_full_sol;
    std::map<const Elem*,std::vector<double> > OverAll;
    std::map<const Elem*,std::vector<std::vector < double> > > Intensity;
    std::map<const Elem*,std::vector< std::vector<double> > >  sol_eignvectors;
    std::vector<double> es;

  
 protected:

    /*! \copydoc SimulationInterface::convert_variable_name_to_id() */
  virtual ID convert_variable_name_to_id(const std::string& variable_name) const{} ;
 
  //! \copydoc  SimulationInterface::do_init() 
  virtual void 	do_init (void);
  
  //!Do the solve
  virtual void do_solve (void);
 
  //!Parse the options
  virtual void 	parse_options(void);
 

};
 
#endif
