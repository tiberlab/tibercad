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
#include "tensor.h"
#include "ModelOptions.h"
#include "AtomisticStructure.h"
#include "TypeDefs.h"
#include "Material.h"

class Material;


//!An internal atom class for Atomistic Generator. 
/*! It differs from normal Atom class as it contains 
 * a flag member useful for some internal routines
 */
class AtomAG
{
public:
  AtomAG(); 
  std::string specie;
  Tensor1 position;
  unsigned int flag; //Default is 0
  
};



class AtomisticGenerator 
{

public:

  AtomisticGenerator(AtomisticStructure* const as);

  ~AtomisticGenerator(void);

  static AtomisticGenerator* create(AtomisticStructure* const as);

   //Tolerance defined internally for casting and comparison
  static const double tol;

  //! Initialize structure informations
  void do_init(void);


  //! Lattice constants  
  double ax, ay, az;

  //! Set primitive vectors, depending on lattice name
  void set_lattice_type(const std::string lattice_name);

 
private:

  //! Primitive vectors in real space, stored by columns in a 3x3 matrix
  Tensor2Gen _prim_vec;

//! 

  //! Name of lattice type
  std::string _lattice_type;

 //! Name of basis type
  std::string _basis_type;

   //Atomic basis arrays
  std::vector<AtomAG> _crystal_basis;

  //! AtomisticStructure instance which invoked AtomisticGenerator 
  //Note: pointer is constant, variable pointed is not constant
  AtomisticStructure*  _as;

  //! Set the atomic basis for the lattice
  void set_crystal_basis(const std::string basis_name, const std::string specie1 = "not_specified", const std::string specie2 = "not_specified", double u = 0.0);

//! Build miller indexes in primitive basis (In some lattices like ZincBlende Miller indexes are not defined
//! on Wiegner-Seitz Cell). Cut_planes are 3X1 Miller indexes arrays stored by columns
void set_prim_miller(Tensor2Gen cut_planes);

//! Cut planes as reciprocal space indexes in reciprocal primitive basis
Tensor2Gen _prim_miller;



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







