using namespace std;

#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>

// Basic include files needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "gmv_io.h"
#include "linear_implicit_system.h"
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

#include "stiffness.h"
#include "rotated_crystal.h"
#include "tensor.h"
#include "GMVIO_cell.h"

#include "piezoelectricity.h"

#include "mesh_data.h"

#include "tecplot_io.h"
#include "tecplot_IO_cell.h"
//------------------------------------------------------------------------------







//! A class to calculate strain and shape
class Macrostrain
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






   //data structure type


  struct options {


    //------------------------------------

    //----------------- numerical options -------------------------------------
    string                mesh_input_file;

    unsigned int          max_r_steps;       
    int                   uniform_refinement; 
     
    double                refine_fraction ; 
    double                coarsen_fraction ;  
    unsigned int          max_ref_level  ;   
    double                tolerance     ;    
  
    unsigned int          max_shape_steps;   
    //--------------------------------
    bool                  grown_on_substrate;
  
    bool                  periodicity[3]; 
    unsigned int          substr_mat;
    //-------------------------------
    
    vector<double> fixed_point1 ; //x,y,z
    vector<double> fixed_point2 ;
    vector<double> fixed_point3 ;
    
 

    unsigned int           number_of_regions;

    //-------------------------------------------------------------------------
    //atomic displacements
  
    bool calculate_atom_displacements;
    
    string   atom_structure_filename;
    
    string   atom_displacements_filename;  
    //-----------------------------------------------------------------------
    //output options
    
    bool intermediate_output;
  
    
    string  output_type;

    //------------------------------------------------------------------------  
  };

  struct  add_variable 
  {
    string name ;
    Elem * element;
    unsigned int dof_number;
    bool lat_cons;
    unsigned int index1;
    unsigned int index2;
    
  } ;
  

  
 
  //---------------------------------------------------------------------
  
  //! Constructor
  /*!
    \param opt options
    \mesh  reference to a mesh object
  */
  Macrostrain(const options& opt,  Mesh& mesh); //gets mesh 
  

  //----------------------------------------------------------------------

  //!set information about crystal orientation and elastisity moduli
  /*!
    \param C_tensor_in vector of objects that can provide elasticity moduli information
    \param crystal_in ector of objects that can provide information about crystals
  */

  void define_strain_parameters(const std::vector<stiffness>&        C_tensor_in,
				const std::vector<rotated_crystal>&  crystal_in);

  

  //--------------------------------------------------------------------
  //! passes a number of substrate boundary condition
  void define_substrate_bc(unsigned int substrate_bc_number);

  //--------------------------------------------------------------------
  
  void define_piezo_moduli(std::vector<Piezoelectricity>&  piezo_in);
  //---------------------------------------------------------------------

  //! passes a reference to a boundary conditions map  
  void define_BC_map (const map <unsigned int , vector<unsigned int> > & bc_cond  );

  //---------------------------------------------------------------------
  //! passes a value for the stress 
  /*!
    \param stress_map is a map between boundary condition number and stress value
  */
  void define_stress_value (const map <unsigned int, double> & stress_map_in);


  void assign_mesh_data(MeshData& mesh_data_in);
  

  //---------------------------------------------------------------------
  /*
    Static function that assembles the linear system matrix for Libmesh
  
   */
  static void assemble_strain_matrix(EquationSystems& es,
				     const std::string& system_name);
  //---------------------------------------------------------------------

  /*
    method that solves strain problem 
  */
  void solve();
  //---------------------------------------------------------------------


  Tensor2Sym get_strain(const Elem* el, bool crystal_system = false); //calculate strain
  
  Tensor1 get_piezopolarization(const Elem* el);
  //---------------------------------------------------------------------
  void output_strain(std::string filename ); //output strain for gmv

  void output_piezo(std :: string filename); //output piezo for gmv


  void output_add_strain_variables(string filename); // output lattice matching parameters
  //---------------------------------------------------------------------

  void output_materials(std :: string filename); //output materials

  //----------------------------------------------------------
 
  Mesh* get_mesh(); //get pointer to the mesh
  


  ~Macrostrain();


  



 private:

  


  MeshData*  meshdata;


  EquationSystems*   equation_systems; //pointer to the equation system

  
  //---------------------------------------------------------------------
  /*
    STATIC Pointers to the objects that are necesary to assemble strain problem matrix

  */
  static std:: vector<stiffness>*        C_tensor_temp;
  static std:: vector<rotated_crystal>*  crystal_temp;
  static std:: vector <int>*             material_of_elem_temp;
  static std:: vector <Tensor2Sym>*      eps0_of_elem_temp;

  static std:: vector<add_variable>*     add_var_temp; 
   
  static std::vector<unsigned int>*      zero_set_dofs_temp;

  //---------------------------------------------------------------------------------------------------
  /*
    DYNAMICAL objects that are necesary to assemble strain problemmatrix
  */
  std::vector<stiffness>        C_tensor; //Elasticity tensor of a material     
  std::vector<rotated_crystal>  crystal;  //crystalographic information of a material 

  std::vector<int>              material_of_elem; // material of an element
  std::vector<Tensor2Sym>       eps0_of_elem;     // eps of an element from previous iteration


  std::vector<unsigned int>     zero_set_dofs;    //DOFs number that have to be set to zero
  //----------------------------------------------------------------------------------------------------

  std::vector< std:: vector <double> > u_node; //displacements of nodes with respect to the non-deformed mesh

  void init_u_node(); //initializes list of node displacements u = 0.0

  void update_u_node(); //updates list of node displacements

  //---------------------------------------------------------------------

  std::vector<Piezoelectricity> piezo; //piezoelectricity constants of a meterial   

  std :: map <const Elem*, unsigned int > elem_numbers;  //map between element pointers and their numbers

  std :: map <const Node*,  unsigned int > active_node_number;    //map between active node pointers and their numbers

  std :: vector< std :: vector <const Node*> >  nodes_periodic; //dim node list's: each contains list of nodes that periodic b.c
                                                                //must be applied to

  static  std::string uname_vec[3];  // "ux", "uy", uz"
  
  static bool belongs_to_substrate(unsigned int n, const Elem* elem ); //if vertex n of elem belongs to substrate
 
  static inline int delta (int i, int j); //Kronecker detla

  static unsigned int dim; //problem dimension

  static bool grown_on_substrate; //if there is a substrate



  void refer_objects(); //refer static pointers to dinamical objects

  void assemble_material_list(); // create material_of_elem at the beginning


  void update_eps0_list();  // update eps0_of_elem // !check!!!!

  void initialize_eps0_list(); //initialize eps0_of_elem !check!!!!

  void initialize_el_number_map(); //initialize   elem_numbers; 

  void make_nodes_periodic(); //create nodes_periodic

  void apply_periodic_bc(); //create DOF constraints for periodic b.c.

  void apply_antirotation_constraints(); //create DOF constraints that do not allow to a freestanding system to rotate.


  Tensor2Sym calculate_eps_lat_matching(unsigned int material); //calculate latiice matching tensor considering latice constants and strain
  

  
  void move_nodes(); //create new mesh by moving nodes slightly


  
  bool periodicity[3]; 


  static unsigned int substr_mat; //substrate material number
  

  

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
  

  string mesh_input_file;

  bool read_regions_from_mesh;
  //------------------------------------------------------------------------
  //Additional variables----------------------------------------------------
  unsigned int             number_of_add_var;
  
  static  unsigned int             number_of_add_var_static;

  vector<add_variable>  add_var;

  void define_additional_variables();//understand which additional varables are necessary.

  void set_up_additional_dofs();

  static std::vector<unsigned int> add_dofs_vector;

  
  static Tensor2Sym  eps0_var_log;  
  /*
    = 0, this is a fixed entry 
    = 1, this is variable entry
  */
  //-------------------------------------------------------------------

  unsigned int fixed_node1;
  unsigned int fixed_node2;
  unsigned int fixed_node3;
  

  static unsigned int fixed_node1_temp;
  static unsigned int fixed_node2_temp;
  static unsigned int fixed_node3_temp;

  Point fixed_point1 ; //x,y,z
  Point fixed_point2 ;
  Point fixed_point3 ;

  //------------------------------------------------------------------
  //!number of boundary condition that defines substrate
  unsigned int substrate_bc_number;

  static double substrate_lat_const[3];

  static Tensor2Sym substrate_shear;
  
  void update_substrate(); //updates substrate_lattice and substrate_shear

  void init_substrate(); //sets substrate_lattice[3] as a reference material

  unsigned int get_number_of_the_fixed_node(Point point); //returns node number, closest to point 

  void Macrostrain::create_substate_nodes_set();
 
  void Macrostrain::update_substrate_nodes_set();
  //------------atomic description------------------

  std::vector<atom> atom_structure;  

  void read_atom_structure(const std::string filename);
  


  void write_atom_displacements(const std::string filename);

  bool calculate_atom_displacements;
 
  string atom_structure_filename;
 
  string atom_displacements_filename;  

  string atom_output_type;

  bool output_strain_on_atoms;

  //-------------------------------------------------



  bool may_belong_to_element(const Elem* element, Point& point);

  unsigned int find_nearest_node(Point& point);


  void define_fixed_nodes();


  //-------------------------------------------------------------------

  bool intermediate_output; // to perform output after each mesh refinement and volume relaxation step

  //-------------------------------------------------------------------


  map<const Elem*, set<unsigned int> > substrate_faces; // map between Elem* and set of substrate faces

  map <unsigned int , vector<unsigned int> >   boundary_cond_nodes; //map between b.c. number and a set of nodes 


  map <const Elem*, map <unsigned int, double>  >   boundary_cond_elem; //map between Elem  and (map between side and stress value)

  static map <const Elem*, map <unsigned int, double>  >* boundary_cond_elem_temp; //static pointer to boundary_cond_elem


  map <unsigned int, double>   stress_values; //map between stress number and stress values
  

  set <unsigned int> substrate_nodes; //contains nodes that belong to substrate
  static set <unsigned int>* substrate_nodes_temp; //static pointer to substate_nodes
   

  void create_bondary_conditions_map(); // create map boundary_cond_elem from boundary_cond_nodes;
   

  void update_bondary_conditions_map(); // update boundary_cond_elem due to mesh refinement ;

  
  string output_type; //"GMV" or "tecplot"


  //-------------------------------------------------------------------
 

};


