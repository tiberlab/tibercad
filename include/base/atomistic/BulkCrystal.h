#ifndef _BULKCRYSTAL_H_
#define _BULKCRYSTAL_H_

#include "AtomisticBasis.h"
#include "ModelOptions.h"
#include "Material.h"
#include "tensor.h"

class BulkCrystal: public AtomisticBasis
{
  public:

  static BulkCrystal* create(const Material* mat, const ModelOptions& options = ModelOptions());

  void init(void);

  /*! Return the vector of Atoms defining the crystal basis
   * unrotated in xyz basis
   */
  const std::vector<Atom>& get_basis(void);

  /*! Return the vector of Atoms defining the crystal basis
   * Already rotated in xyz basis
   */
  const std::vector<Atom>& get_rotated_basis(void);

  /*! Return the vector of Atoms defining the crystal basis
   *  in primitive vectors units
   */
  const std::vector<Atom>& get_lattice_basis(void);

  //! Return string specifying lattice type (ex: 'fcc', 'hexagonal')
  const std::string get_lattice_type(void);

  //! Return the lattice constant as a vector [a1,a2,c]
  const std::vector<double>& get_lattice_constant(void);

  //! Return primitive vectors (3x3 tensor, stored columnwise)
  const Tensor2Gen& get_prim_vec(void);

  //! Return primitive vectors (3x3 tensor, stored columnwise)
  const Tensor2Gen& get_rotated_prim_vec(void);

  //! Return rotation matrix
  const Tensor2Gen& get_rotation(void);

  protected:
    
  private:

  /*! Constructor accepts a material and additional options
   * Option from material and additional optional are merged 
   * internally
   */
  BulkCrystal(const Material* mat, const ModelOptions& options);

  //! Contains the options, as given when create is invoked
  ModelOptions _options;

  //!Read useful information from material database
  void read_database();

  //! Reference 
  const Material* _mat;

  //lattice constants
  std::vector<double> _lattice_constant;
 
  //angles
  std::vector<double> _angles;
 
  //! Atomic basis arrays (non rotated)
  std::vector<Atom> _basis;

  //! Atomic basis arrays 
  std::vector<Atom> _rotated_basis;

  //! Atomic basis in primitive vectors units
  std::vector<Atom> _lattice_basis;

  //! Name of lattice type
  std::string _lattice_type;

  //! Set primitive vectors, according to lattice
  void set_prim_vec(void);

  //! Set rotation matrix and apply rotation where needed
  void build_rotation(void);

  //! Primitive vectors (lattice constant included). Columnwise
  Tensor2Gen _prim_vec;

  //! Rotated primitive vectors
  Tensor2Gen _rotated_prim_vec;

  //!Rotation matrix
  Tensor2Gen _rotation;

};


inline
const std::vector<Atom>&
BulkCrystal::get_basis(void)
{
  return _basis;
}


inline
const std::vector<Atom>&
BulkCrystal::get_lattice_basis(void)
{
  return _lattice_basis;
}


inline
const std::vector<Atom>&
BulkCrystal::get_rotated_basis(void)
{
  return _rotated_basis;
}


inline
const std::string
BulkCrystal::get_lattice_type(void)
{
  return _lattice_type;
}


inline
const std::vector<double>&
BulkCrystal::get_lattice_constant(void)
{
  return _lattice_constant;
}


inline
const Tensor2Gen& 
BulkCrystal::get_prim_vec(void)
{
  return _prim_vec;
}

inline
const Tensor2Gen& 
BulkCrystal::get_rotated_prim_vec(void)
{
  return _rotated_prim_vec;
}


inline
const Tensor2Gen&
BulkCrystal::get_rotation(void)
{
  return _rotation;
}

#endif // _BULKCRYSTAL_H_
