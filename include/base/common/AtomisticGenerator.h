#ifndef _ATOMISTICGENERATOR_H_
#define _ATOMISTICGENERATOR_H_

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <map>
#include <set>
#include "AtomisticStructure.h"
#include "Atom.h"
#include "tensor.h"
#include "ModelOptions.h"
#include "TypeDefs.h"
#include "Material.h"
#include "mesh.h"
#include "mesh_base.h"
#include "elem.h"
#include "point.h"
#include "node.h"

//class Material;
	
class AtomisticGenerator 
{

public:

  AtomisticGenerator(AtomisticStructure* const as);

  ~AtomisticGenerator(void);

  static AtomisticGenerator* create(AtomisticStructure* const as);

  //! Tolerance defined internally for casting and comparison
  static const double tol;

  //! Scaling value respect to TiberCAD units (usually Amstrong instead of micron)
  static const double scale;

  //! Initialize structure informations
  void do_init(void);


  //! Lattice constants  
  double ax, ay, az;

  //! Set primitive vectors, depending on lattice name
  void set_lattice_type(const std::string lattice_name);

  //Print atom_basis in xyz file (for debugging)
  void print_basis(std::vector<Atom> &basis, const std::string filename);
 
private:

  //! Change atom species according to regions
  void change_specie();

  //! Primitive vectors in real space, stored by columns in a 3x3 matrix
  Tensor2Gen _prim_vec;

  //! Conventional cell vectors in absolute basis
  Tensor2Gen _conv_vect;
  //! Conventional cell vectors in primitive real basis
  Tensor2Gen _conv_prim;

  //! Conventional cell basis atoms
  std::vector<Tensor1> conv_lattice_basis;

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
  
  //! Conventional cell lattice points
  std::vector<Tensor1> _conv_lattice_basis;
  
  //! Supercell atom basis points
  std::vector<Atom> _super_basis;

  //! Final structure basis
  std::vector<Atom> _structure_basis;

  //!Supercell lattice points
  std::vector<Tensor1> _super_lattice;
  
  //!Supercell periodicity vectors
  Tensor2Gen _period;

  //! AtomisticStructure instance which invoked AtomisticGenerator 
  //Note: pointer is constant, variable pointed is not constant
  AtomisticStructure*  _as;

  //! List of elements covered by structure
  std::vector<Elem*> _structure_elements;

  //! Set the atomic basis for the lattice
  void set_crystal_basis(const std::string basis_name, const std::string specie1 = "not_specified", const std::string specie2 = "not_specified", double u = 0.0);

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
  void make_supercell(double l1, double l2, double l3, bool preserve_basis, bool preserve_conv);

  //! Function with common operations for building atomistic structure (local calls to build1D, build2D, build3D)
  void build();

  //! Function which builds 1D atomistic structure (one conventional cella in 0 directions)
  void build1D();


  //Some data manipulation function useful only in this class
  
  //Calculate a reciprocal basis from a real basis
  Tensor2Gen reciprocal(Tensor2Gen real_basis);
  
  //Find greater common denominator between two integers
  int gcd(int a, int b);

  //Reduce vector, dividing its members for GCD (cast from double to int temporary)
  Tensor1 reduce_vector(Tensor1 v);

  //Reduce vectors contained in a 2D tensor
  Tensor2Gen reduce_vector(Tensor2Gen a);

  //Comparison with tolerance
  int compare_tol(double a, double b);

  //Casting from double to int, checking that fractional part is minor than an internal  tolerance
  int double_to_int_cast_checked(double a);

  //Reduce double to nearest integer and checked whatever difference is greater than internal tolerance
  double double_to_int_value_checked(double a);

  //same on entire tensor (in place)
  void double_to_int_value_checked(Tensor1& a);

  //Scale a vector with fractional parts to all integer vector
  void scale_to_int(Tensor1& a); 

  //Same work with 3 vectors composing a tensor2Gen
  void scale_to_int(Tensor2Gen& a);

};






#endif // _ATOMISTICSTRUCTURE_H_







