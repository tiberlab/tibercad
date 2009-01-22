#ifndef _ATOMISTICGENERATOR_H_
#define _ATOMISTICGENERATOR_H_

//--------------------------------------------------------------------------------------------

#include "Atom.h"
#include "tensor.h"
#include "ModelOptions.h"
#include "TypeDefs.h"
#include "Material.h"
#include "Database.h"
#include "Alloy.h"
#include "mesh.h"


//! A class for building Atomistic Structure from mesh informations
/*!
 *Atomistic Generator can create 1D, 2D and 3D atomistic structure
 *(directions confinement define structure periodicity), according
 *to mesh informations. Material parameters (kind of lattice, atomic basis, ecc.)
 *are read from material files
 */


//forward declaration
class AtomisticStructure;
class BondMap;




class AtomisticGenerator
{

public:

  AtomisticGenerator(void);

  virtual   ~AtomisticGenerator(void);

  //! Initialize structure informations
  void do_init();

  static AtomisticGenerator* create(AtomisticStructure* const as, unsigned int dimension);

  //! Tolerance defined internally for casting and comparison
  static const double tol;

  //! Scaling value respect to TiberCAD units (usually Amstrong instead of micron)
  double scale;

  //! Set primitive vectors, depending on lattice name
  void set_lattice_type(const std::string lattice_name);

  //Print atom_basis in xyz file (for debugging)
  void print_basis(std::vector<Atom> &basis, const std::string filename);


protected:

  //lattice constants
  double _lattice_constant[3];

  //! Change atom species according to regions
  void cut_and_change_specie(std::string preserve);

  //  //! Fast bond map generation, suitable for both
  //   unsigned int** fast_bond_map(const std::vector<Atom> &basis,
  //		   Tensor1& edge_min, Tensor1& edge_max, Tensor2Gen& period);

  //! Primitive vectors in real space, stored by columns in a 3x3 matrix
  Tensor2Gen _prim_vec;

  //! Conventional cell vectors in absolute basis
  Tensor2Gen _conv_vect;
  //! Conventional cell vectors in primitive real basis
  Tensor2Gen _conv_prim;

  //! Conventional cell basis atoms
  std::vector<Tensor1> _conv_lattice_basis;

  //!Rotation tensor
  Tensor2Gen _rotation;

  //! Localo origin
  Tensor1 _local_origin;

  //! Name of lattice type
  std::string _lattice_type;

  //! Name of basis type
  std::string _basis_type;

  //! Atomic basis arrays
  std::vector<Atom> _crystal_basis;

  //! Supercell atom basis points
  std::vector<Atom> _super_basis;

  //! Final structure basis
  std::vector<Atom> _structure_basis;

  //!Supercell lattice points
  std::vector<Tensor1> _super_lattice;

  //! Supercell conventional cell edges points
  std::vector<Tensor1> _super_conv;

  //! Dimensionality of the system (1, 2 or 3)
  unsigned int _dim;

  //! Lenght of supercell in conventional cells units
  int _conv_cells_supercell_lenght[3];

  //! An internal instance of BondMap, to make passivation and final bond map
  BondMap* _bondmapobject;

  //! Reference material
    Material* _reference_material;

    //! Database of reference material
    Database _reference_material_db;

  //! Missing super_conv (vector of conventional cells edges).
  //!If it will be needed remember to uncomment proper lines in make_supercell!!!!!!!!!!!!!!!!


  //Bond map generation
  void bond_map_gen(std::vector<Atom> &basis);



  // Build cutoff distancies map
  void set_cutoff();


  //Passivation routine for bulk structures (periodicization is achieved in derived classes, in hydrogenation routine)
  void passivate_cluster(std::vector<Atom> &basis);


  //! Real passivation routine (implemented in derived classes, it takes in account periodicity)
  virtual void passivate() = 0;


  //!Supercell periodicity vectors
  Tensor2Gen _period;

  //Map for cutoff parameters
  std::map<std::string, double> _cutoff;

  //! AtomisticStructure instance which invoked AtomisticGenerator
  //Note: pointer is constant, variable pointed is not constant
  AtomisticStructure*  _as;

  //! List of elements covered by structure
  std::vector<Elem*> _structure_elements;

  //! Set the atomic basis for the lattice (This function is no longer used!)
  //void set_crystal_basis(const std::string basis_name, const std::string specie1 = "not_specified", const std::string specie2 = "not_specified", double u = 0.0);

  //! Build miller indexes in primitive basis (In some lattices like ZincBlende Miller indexes are not defined
  //! on Wiegner-Seitz Cell). Cut_planes are 3X1 Miller indexes arrays stored by columns
  void set_prim_miller(Tensor2Gen cut_planes);

  //! Cut planes as reciprocal space indexes in reciprocal primitive basis
  Tensor2Gen _prim_miller;

  //! Setting growth conventional cell vectors (in primitive vectors basis)
  void make_conv_cell();

  //! Filling conventional growth cell with basis atoms
  void make_conv_basis();

  //! A function to build supercells
  void make_supercell(double l1, double l2, double l3);

  //! Virtual function for building up the structure.
  virtual void build() = 0;

  //! Parsing of atomistic infos to build lattice and basis vectors
  void parse_parameters(const Material* mat);


  //Some data manipulation function useful only in this class




  //Calculate a reciprocal basis from a real basis
  static Tensor2Gen reciprocal(Tensor2Gen real_basis);

  //Find greater common denominator between two integers
  static int gcd(int a, int b);

  //Reduce vector, dividing its members for GCD (cast from double to int temporary)
  static Tensor1 reduce_vector(Tensor1 v);

  //Reduce vectors contained in a 2D tensor
  static Tensor2Gen reduce_vector(Tensor2Gen a);

  //Comparison with tolerance
  static int compare_tol(double a, double b);

  //Casting from double to int, checking that fractional part is minor than an internal  tolerance
  static int double_to_int_cast_checked(double a);

  //Reduce double to nearest integer and checked whatever difference is greater than internal tolerance
  static double double_to_int_value_checked(double a);

  //same on entire tensor (in place)
  static void double_to_int_value_checked(Tensor1& a);

  //Scale a vector with fractional parts to all integer vector
  static void scale_to_int(Tensor1& a);

  //Same work with 3 vectors composing a tensor2Gen
  static void scale_to_int(Tensor2Gen& a);

};






#endif // _ATOMISTICGENERATOR_H_







