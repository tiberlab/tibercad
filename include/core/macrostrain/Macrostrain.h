#ifndef _MACROSTRAIN_H_
#define _MACROSTRAIN_H_


#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

// Basic include files needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "gmv_io.h"
#include "equation_systems.h"

#include "getpot.h"

// For mesh refinement
#include "mesh_refinement.h"
#include "error_vector.h"
#include "kelly_error_estimator.h"

// Define the Finite Element object.
#include "fe.h"
#include "elem.h"
// Define Gauss quadrature rules.
#include "quadrature_gauss.h" 

// Define useful datatypes for finite element
// matrix and vector components.
#include "sparse_matrix.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "dof_map.h"

#include "fe_interface.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>



#include <fstream>
#include <iomanip>

//------------------------------------------------------------------------------

#include "Stiffness.h"
#include "RotatedCrystal.h"

#include "GMVIO_cell.h"

#include "Piezoelectricity.h"

#include "mesh_data.h"

#include "tecplot_io.h"
#include "tecplot_IO_cell.h"
//------------------------------------------------------------------------------

#include "SimulationInterface.h"
#include "StrainSimulation.h"
#include "Device.h"
#include "MacrostrainModel.h"

class TiberLinearSystem;

//! A class to calculate strain and shape
class Macrostrain : public StrainSimulation
/*!
  Class uses the macroscopic strain theory
*/

{
 public:





  struct atom
  {
    int mat_number;
    int type;
    Elem* element;
    Point relative_point;

  };








  struct  add_variable 
  {
    std::string name ;
    Elem * element;
    unsigned int dof_number;
    bool lat_cons;
    unsigned int index1;
    unsigned int index2;
    
  } ;
  
  //!structure that contains both crystalographic information of a material and Elasticity tensor of a material
  struct strain_param
  {
    Stiffness* C_tensor;
    RotatedCrystal* crystal;
  };
 
  
  //!Constructor
   Macrostrain(void );


  

  //---------------------------------------------------------------------
  /*
    Static function that assembles the linear system matrix for Libmesh
  
   */
  static void assemble_strain_matrix(EquationSystems& es,
				     const std::string& system_name);
  //--------------------------------------------------------------------
  

  //!get in crystal strain system
  /*!
    \param el pointer to the element
    \param quadratur_point quadratur point that belongs to the element
  */
  Tensor2Sym get_strain_crystal(const Elem* el, const Point& quadratur_point ); 
  
  //!get in crystal strain system
  /*!
    \param el pointer to the element
  */
  Tensor2Sym get_strain_crystal(const Elem* el);
  
 //Tensor2Sym get_stress_crystal(const Elem* el);
 

  //!get polarization (piezo)
  /*!
    \param el pointer to the element
  */
  Tensor1 get_piezopolarization(const Elem* el);


  //!get built-in polarization (piezo)  
 
  /*!
    \param el pointer to the element
    \param quadratur_point quadratur point that belongs to the element
  */
  Tensor1 get_built_in_polarization(const Elem* el, const Point& quadratur_point ); 

 
 
  //!output piezo for gmv
  void output_piezo(std :: string filename); 

  //! output lattice matching parameters
  void output_add_strain_variables(std::string filename); 
  //---------------------------------------------------------------------



  //----------------------------------------------------------
 
  Mesh* get_mesh(); //get pointer to the mesh
  


  
  


  virtual ~Macrostrain();


   
  static Macrostrain* create(void);
      

 


  /*! \copydoc SimulationInterface::create_boundary_model() */
  virtual BoundaryProperties*
    create_boundary_model(const ModelOptions& options) const
    throw (ModelErrorException);
  


 
     
     

 
 

  //!will be removed in future
  void write_atom_potential();

  //! Preliminary check to see if a point could belong to an element (much faster than exact calculation)
  static bool may_belong_to_element(const Elem* element, Point& point);

 private:


  //!if true then all the fixed points are zero displacement
  bool fix_all_fixed_points;


  //! reallocate system matrix
  void reallocate_matrix(void);

  bool _is_reallocated;
  bool _preallocate;

  //!true the constrains have to be applied
  bool _first_run;

  //!pointer to a drift-diffusion object that is used to get electric-field  data 
  SimulationInterface* poisson_equation;

  //!pointer to the equation systems 
  EquationSystems*   equation_systems; 

  

  //!pointer to the  system used in the simulation 
  TiberLinearSystem* my_system;
  
 

  //!name of my system
  std::string system_name;

  


  //!pointer to the mesh
  Mesh*  mesh;

  //!calculate strain
  Tensor2Sym get_strain(const Elem* el, bool crystal_system = false); 
 
  //! calculate the result_strain map
  void calculate_result_elem_strain_map();

  
  
  //---------------------------------------------------------------------
 
  //! static pointer to this object
  static Macrostrain* static_this; 



  //! Substrate material crystal
  const RotatedCrystal* substrate_crystal;

 
  //! strain tensor of an element from previous iteration (calculation system)
  std::vector<Tensor2Sym>       eps0_of_elem;    


 
 
  //!displacements of nodes with respect to the non-deformed mesh
  std::map< const Node*, std:: vector <double> > u_node; 

  //!initializes list of node displacements u = 0.0
  void init_u_node(); 

  //!updates list of node displacements
  void update_u_node(); //updates list of node displacements

  
 
  //!map between element pointers and their numbers
  std :: map <const Elem*, unsigned int > elem_numbers; 
  
 
  //!dim node list's: each contains list of nodes that periodic b.c. must be applied to
  std :: vector< std :: vector <const Node*> >  nodes_periodic; 

  //! {"ux", "uy", uz"}
  std::string uname_vec[3];  


  //!if vertex n of elem belongs to substrate
  bool belongs_to_substrate(unsigned int n, const Elem* elem ); //if vertex n of elem belongs to substrate

  //!Kronecker detla
  inline int delta (int i, int j); //Kronecker detla

  //!problem dimension
  unsigned int dim; //problem dimension

  //!if there is a substrate
  bool grown_on_substrate; 

  //!  name of the substrate boundary condition;
  std::string substrate_name;

  //!refer static pointers to dinamical objects
  void refer_objects(); //refer static pointers to dinamical objects



  //! update eps0_of_elem 
  void update_eps0_list();  // update eps0_of_elem // !check!!!!

  //!initialize eps0_of_elem 
  void initialize_eps0_list(); 

  //!initialize   elem_numbers;
  void initialize_el_number_map(); 

  //!create nodes_periodic
  void make_nodes_periodic(); 


  //!create DOF constraints for periodic b.c.
  void apply_periodic_bc(); 

  //!create DOF constraints that do not allow to a freestanding system to rotate.
  void apply_antirotation_constraints(); 

  //!calculate latiice matching tensor considering latice constants and strain
  Tensor2Sym calculate_eps_lat_matching(const Elem* elem); 
  

  //!create new mesh by moving nodes slightly
  void move_nodes();

 
  //!writes in a file all strain components (only for debug purposes)
  /*!
    \filename name of the file for output
  */
  void output_strain(std::string filename );
  
  bool periodicity[3]; 


  

  //----device boundaries---------------
  double min_coord[3];
  double max_coord[3];
  //------------------------------------

  //----------------- numerical options -------------------------------------
  unsigned int          max_r_steps;       
  int                   uniform_refinement; 
     
  double                refine_fraction ; 
  double                coarsen_fraction ;  
  unsigned int          max_ref_level  ;   
  double                tolerance     ;    
  
  unsigned int          max_shape_steps;   
  //------------------------------------------------------------------------
  // mesh reading
  

  std::string mesh_input_file;

  bool read_regions_from_mesh;
  //------------------------------------------------------------------------
  //Additional variables----------------------------------------------------
  unsigned int             number_of_add_var;
  
 

  std::vector<add_variable>  add_var;

  void define_additional_variables();//understand which additional varables are necessary.

  void set_up_additional_dofs();

  std::vector<unsigned int> add_dofs_vector;

  
   Tensor2Sym  eps0_var_log;  
  /*
    = 0, this is a fixed entry 
    = 1, this is variable entry
  */
  //-------------------------------------------------------------------

  unsigned int fixed_node1;
  unsigned int fixed_node2;
  unsigned int fixed_node3;
 
 

  Point fixed_point1 ; //x,y,z
  Point fixed_point2 ;
  Point fixed_point3 ;
 

  //------------------------------------------------------------------
  //!number of boundary condition that defines substrate
  unsigned int substrate_bc_number;

  double substrate_lat_const[3];

  Tensor2Sym substrate_shear;

  //!updates substrate_lattice and substrate_shear
  void update_substrate();

  //!sets substrate_lattice[3] as a reference material
  void init_substrate();

  //!returns node number, closest to point 
  unsigned int get_number_of_the_fixed_node(Point point); 


 
  //------------atomic description------------------

  std::vector<atom> atom_structure;  

  void read_atom_structure(const std::string filename);
  


  void write_atom_displacements(const std::string filename);

 

  bool calculate_atom_displacements;
 
  std::string atom_structure_filename;
 
  std::string atom_displacements_filename;

  std::string atom_potential_filename;
  

  std::string atom_output_type;

  bool output_strain_on_atoms;

  //-------------------------------------------------



 

  unsigned int find_nearest_node(Point& point);


  void define_fixed_nodes();


  //-------------------------------------------------------------------
  //! to perform output after each mesh refinement and volume relaxation step
  bool intermediate_output;
  //-------------------------------------------------------------------


  //!node id's that belong to substrate
  std::set <const Node*> substrate_points;

 

  //!"GMV" or "tecplot"
  std::string output_type; 

  //!
  bool element_on_boundary(const Elem* element);
 

  //! non-static method that actually does matrix assembling 
  void do_assemble(EquationSystems& es,  const std::string& system_name);



  //! calculates \f$ \mathop{\rm{max}}_{\alpha, n}|u_{\alpha}^n - v_{\alpha}^n|  \f$
  double norm_of_difference(NumericVector<Number>& solution1, NumericVector<Number>& solution2);

  //! Preapare all 6 components of the strain tensor for output
  void prepare_strain_data_for_output( std::vector<std::string>& eps_names, std::vector<double>& eps_data );

  //! Preapare all 6 components of the stress tensor for output
  void prepare_stress_data_for_output( std::vector<std::string>& stress_names, std::vector<double>& stress_data );

  //! Preapare all 3 components of the polarization vector for output
  void prepare_polarization_data_for_output( std::vector<std::string>& polariz_names, std::vector<double>& polariz_data );


  //! recalculate partial the derivatives of the basis functions for the zero-gradient boundary condition
  /*!
    \param deriv_vectors input: initial derivatives; output: recalculated derivatives
    \param normal normal vector
  */
  void adjust_derivatives(Tensor1& deriv_vectors, const Point& normal); 
 
 protected:

 

  virtual void do_init(void);

  virtual void do_solve(void);
  
  virtual void parse_options(void);


  

  
  /*! 
    \copydoc SimulationInterface::build_elemental_results()
    The variables are: "strain", "polarization"
    This means strain tensor components:
    \f$ \varepsilon_{xx}, \varepsilon_{yy},\varepsilon_{zz}, \varepsilon_{xy}, \varepsilon_{xz}, \varepsilon_{yz}\f$,
    and polarization vector \f$ \bf P \f$ components. 
    \f$ x, y,  z\f$ refer to calculation coordinate system.
   */
  virtual void build_elemental_results(const std::set<std::string>& variables,
				       std::vector<double>& results, std::vector<std::string>& legend) ;
   
 


  
 



};


//-------------------------------------------------------------------
inline bool Macrostrain::element_on_boundary(const Elem* element)
{
  bool result = false;

  const Mesh& mesh = equation_systems->get_mesh();
  
    
  unsigned int n_sides ; 

  if ( dim > 1 ) 
    n_sides = element->n_sides();
  else
    n_sides = element->n_nodes();


  for (short i = 0; i < n_sides; i++)
    {
      Elem* el1 = element->neighbor(i);
      if ( (el1 == NULL)  ) 
	  result = true;
      else
	if (!( el1 -> active() ))
	  result = true;
	  
      if (result) break;
	
      
    }

 
  return(result);
  
}


inline
Tensor2Sym
Macrostrain::get_strain_crystal(const Elem* el)
{
  return get_strain_crystal(el, el->centroid());
}

inline
Macrostrain* Macrostrain::create(void) 
{
  return new Macrostrain();
}

#endif
