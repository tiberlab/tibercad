// $Id$

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



// forward declarations
class AtomisticStructure;
class BondMap;
class BulkCrystal;



//! A class for building Atomistic Structure from mesh informations
/*!
 *Atomistic Generator can create 1D, 2D and 3D atomistic structure
 *(directions confinement define structure periodicity), according
 *to mesh informations. Material parameters (kind of lattice, atomic basis, ecc.)
 *are read from material files
 */
class TBDLLOCAL AtomisticGenerator
{

public:

  virtual  ~AtomisticGenerator(void);

  //! Initialize structure informations
  void do_init();

  //! Create an object of AtomisticGenerator for the right dimensionality
  static AtomisticGenerator* create(AtomisticStructure* const as, unsigned int dimension);

  //! Tolerance defined internally for casting and comparison
  static const double tol;

  //! Scaling value respect to TiberCAD units (usually Amstrong instead of micron)
  double scale;

  //! Print atom_basis in xyz file (for debugging)
  void print_basis(std::vector<Atom> &basis, const std::string filename);

  //! Copy back information to AtomisticStructure
  void finalize(void);

  //! associate element to each atom
  /*!
   * This iterates on _super_basis
   */
  void associate_elements(const std::set<ID>& reg_ids);
  
  //! cut the structure (only flags atoms)
  void cut(const std::set<ID>& reg_ids, const std::string preserve = "none");

  void dorestrict(bool passivation = true);
  
  //! assign the correct specie to each atom
  /*!
   * iterates on _super_basis and creates _structure_basis with active atoms only
   */
  void assign_species(void);
  

protected:

  explicit AtomisticGenerator(AtomisticStructure* const as);

  //! Change atom species according to regions
  void cut_and_change_specie(std::string preserve);


  //! Conventional cell vectors in absolute basis
  Tensor2Gen _conv_vect;

  //! Conventional cell vectors in primitive real basis
  Tensor2Gen _conv_prim;

  //! Lattice points in conventional basis
  std::vector<Tensor1> _conv_lattice;

  //! Conventional cell basis atoms
  std::vector<Atom> _conv_basis;

  //! Localo origin
  Tensor1 _local_origin;
  
  //! global translation vector
  Tensor1 _translation;
  
  //! cell translation vector
  Tensor1 _cell_translation;
  
  //! Name of basis type
  std::string _basis_type;

  //! Supercell atom basis points
  std::vector<Atom> _super_basis;

  //! Final structure basis
  std::vector<Atom> _structure_basis;

  //! Dimensionality of the device mesh (1, 2 or 3)
  unsigned int _dim;

  //! Lenght of supercell in conventional cells units
  int _conv_cells_supercell_lenght[3];

  //! An internal instance of BondMap, to make passivation and final bond map
  BondMap* _bondmap;

  //! Reference material
  const Material* _reference_material;

  //! Bond map generation
  void bond_map_gen(const std::vector<Atom>& basis);
  
  //! enlarge supercell vectors for dummy periodicities
  void check_periodic(void);

  //! Actually remove marked atoms from structure
  void remove_atoms(void);

  //! Real passivation routine (implemented in derived classes, it takes in account periodicity)
  void passivate();


  //!Supercell periodicity vectors
  Tensor2Gen _period;

  //! Map for cutoff parameters
  std::map<std::string, double> _cutoff;

  //! AtomisticStructure instance which invoked AtomisticGenerator
  /*!
   * Note: pointer is constant, variable pointed is not constant
   */
  AtomisticStructure*  _as;

  //! Setting growth conventional cell vectors (in primitive vectors basis)
  void make_conv_cell();

  //! Setting the conventional lattice
  void make_conv_lattice();

  void minimal_conv_cell();

  //! Setting the conventional lattice + basis
  void make_conv_basis();

  //! A function to build supercells (basis+lattice filling space)
  void make_supercell(double l1, double l2, double l3);

  //! Function for building up the structure.
  virtual void build(void);

  //! Find the origin of the conventional cell such that all atoms have positive coordinates
  void move_origin(void);


  //! Calculate a reciprocal basis from a real basis
  static Tensor2Gen reciprocal(Tensor2Gen real_basis);

  //! Find greater common denominator between two integers
  static int gcd(int a, int b);

  //! Reduce vector, dividing its members for GCD (cast from double to int temporary)
  static Tensor1 reduce_vector(Tensor1 v);

  //! Reduce vectors contained in a 2D tensor
  static Tensor2Gen reduce_vector(Tensor2Gen a);

  //! Comparison with tolerance
  static int compare_tol(double a, double b);

  //! Casting from double to int, checking that fractional part is minor than an internal  tolerance
  static int double_to_int_cast_checked(double a);

  //! Reduce double to nearest integer and checked whatever difference is greater than internal tolerance
  static double double_to_int_value_checked(double a);

  //! Reduce double to nearest integer and checked whatever difference is greater than internal tolerance
  static void double_to_int_value_checked(Tensor1& a);

  //! Scale a vector with fractional parts to all integer vector
  static void scale_to_int(Tensor1& a);

  //! Scale 3 vectors (as tensor) with fractional parts to all integer vector
  static void scale_to_int(Tensor2Gen& a);

  //! fold atoms into conventional cell
  static bool fold_in_cell(Atom& atom, const Point& orig, const Point& a1, const Point& a2, const Point& a3, bool fold=true);
  //! fold point into conventional cell
  static bool fold_in_cell(Point& p, const Point& orig, const Point& a1, const Point& a2, const Point& a3, bool fold=true);

private:

  AtomisticGenerator(void) {};

  //! Common 0d,1d,2d,3d init operations
  void init_commons();

  //! Build random alloy structure
  void build_random_alloy(void);
  
  //! Calculate environment dependent substitution probability
  double substitution_probability(size_t id, const Specie& sp);

  //! BulkCrystal of the reference material
  BulkCrystal* _bulk;

  //! tells for every atom if it belongs to the structure
  /*!
   * Used for internal purpose
   */
  std::vector<bool> _belong_to_structure;

};





#endif // _ATOMISTICGENERATOR_H_







