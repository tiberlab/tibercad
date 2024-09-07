#ifndef _BULKCRYSTAL_H_
#define _BULKCRYSTAL_H_

#include "AtomisticBasis.h"
#include "ModelOptions.h"
#include "Material.h"
#include "tensor.h"

//! Contains data for a bulk crystal
/*!
 * This class contains all the data of a bulk crystal, including its atomic
 * structure, lattice constants and rotation matrix (from standard crystal
 * directions to simulation coordinates).
 * 
 * Length units are Angstrom.
 *
 */
class BulkCrystal: public AtomisticBasis
{
public:

  //! Creator for BulkCrystal objects
  static BulkCrystal* create(const Material *mat,
                             const ModelOptions &options = ModelOptions());

  //! Initialize the object
  void init(void);

  /*! Return the vector of Atoms defining the crystal basis
   * unrotated in xyz basis
   */
  const std::vector<Atom>& get_basis(void) const;

  /*! Return the vector of Atoms defining the crystal basis
   * Already rotated in xyz basis
   */
  const std::vector<Atom>& get_rotated_basis(void) const;

  /*! Return the vector of Atoms defining the crystal basis
   *  in primitive vectors units
   */
  const std::vector<Atom>& get_lattice_basis(void) const;

  //! Return string specifying lattice type (ex: 'fcc', 'hexagonal')
  const std::string get_lattice_type(void) const;

  //! Return the lattice constants for an orthogonal cell as a vector
  const std::vector<double>& get_ortho_lattice_constants(void) const;

  //! Return primitive vectors in reference coordinates (3x3 tensor, stored columnwise)
  const Tensor2Gen& get_prim_vec(void) const;

  //! Return primitive vectors (3x3 tensor, stored columnwise)
  const Tensor2Gen& get_rotated_prim_vec(void) const;

  //! Return conventional vectors in reference coordinates (3x3 tensor, stored columnwise)
  const Tensor2Gen& get_conv_vec(void) const;

  //! Return conventional vectors (3x3 tensor, stored columnwise)
  const Tensor2Gen& get_rotated_conv_vec(void) const;

  //! Return rotation matrix
  const Tensor2Gen& get_rotation(void) const;

  //! Get the Euler angles
  /*!
   * The Euler angles here are the ones that rotate the crystal from
   * standard orientation to its real orientation in the calculation 
   * coordinate system.
   * We use the convention \f$ Z_\gamma Y_\beta Z_\alpha \f$.
   */
  void get_euler_angles(double& alpha, double& beta, double& gamma) const;

  //! Get lattice matching strain in calculation coordinate system
  /*!
   * This method calculates lattice matching strain, assuming
   * pseudomorphic growth, with respect to a given substrate.
   * we first check if the crystals are compatible. Compatibility
   * is not easy to establish. For now, we check crystal type, and 
   * zb/wz compatibility in terms of growth directions.
   * 
   * \param substrate the substrate BulkCrystal
   * \param strain the lattice matching strain in calculation system
   * \return \c false if the crystals are incompatible  
   */
  bool get_lattice_matching_strain(const BulkCrystal& substrate,
                                   libMesh::RealTensor& strain) const;

  //! Print information on output
  void print_info(void) const;


protected:
private:
  /*! Constructor accepts a material and additional options
   * Option from material and additional optional are merged 
   * internally
   */
  BulkCrystal(const Material* mat, const ModelOptions& options);

  //! Read useful information from material database
  void read_database();

  //! Extract the crystal direction for given Miller indices
  /*!
   * \param dir the direction given as \c "x", \c "y" or \c "z"
   * \param vec the vector in the crystal coordinate system
   * 
   * \return \c true if a valid direction was found
   * 
   * The code will read \c <dir>-growth-direction from options,
   * check whether it is a crystal direction [uvw] or a plane (hkl)
   * or (hkil), and calculate the vector in the crystal system.
  */
  bool extract_crystal_direction(const std::string& dir, Tensor1& vec) const;

  //! Get Miller indices for a given direction
  void get_miller_indices(const Tensor1& vec, std::vector<int>& miller) const;

  //! Get orthogonal vector
  void get_orthogonal_vector(const Tensor1& dir, Tensor1& ortho) const;

  //! Set primitive vectors, according to lattice
  void set_cell_vectors(void);

  //! Set rotation matrix and apply rotation where needed
  /*!
   * This method also recalculates all three directions and
   * the according Miller indices, and the lattice constants
   * along the cartesian axes of the calculation system.
   */
  void build_rotation(void);

  //! Calculate Euler angles for rotation matrix
  void calculate_euler_angles(void);


  //! Contains the options, as given when create is invoked
  ModelOptions _options;

  //! Reference to material
  const Material* _mat;

  //! Lattice constants in Angstrom
  std::vector<double> _lattice_constant {0, 0, 0};
 
  //! Crystal angles, given in degrees
  std::vector<double> _angles {90, 90, 90};
  
  //! Lattice constants for the minimal orthogonal cell, in Angstrom
  /*!
   * These are mainly useful for calculating lattice matching
   * strain, and for matching crystal interfaces.
   */
  std::vector<double> _ortho_lattice_constants {0, 0, 0};
 
  //! Atomic basis arrays (non rotated)
  std::vector<Atom> _basis;

  //! Atomic basis arrays 
  std::vector<Atom> _rotated_basis;

  //! Atomic basis in primitive vectors units
  std::vector<Atom> _lattice_basis;

  //! Name of lattice type
  std::string _lattice_type;

  //! Primitive vectors (lattice constant included). Columnwise
  Tensor2Gen _prim_vec {0};

  //! Rotated primitive vectors
  Tensor2Gen _rotated_prim_vec {0};

  //! Conventional cell vectors. Columnwise
  Tensor2Gen _conv_vec {0};

  //! Rotated conventional cell vectors
  Tensor2Gen _rotated_conv_vec {0};

  //! Rotation matrix
  Tensor2Gen _rotation {1};

  //! The Euler angles \f$\alpha$\f, \f$\beta$\f and \f$\gamma$\f
  /*!
   * The Euler angles here are the ones that rotate the crystal from
   * standard orientation to its real orientation in the calculation
   * coordinate system, given in radians.
   * We use the convention \f$ Z_\gamma Y_\beta Z_\alpha \f$.
   */
  std::vector<double> _euler_angles {0, 0, 0};

  //! Reciprocal lattice vectors, calculated from conventional cell
  std::vector<Tensor1> _reciprocal_lattice = std::vector<Tensor1>(3, Tensor1(0));
};

inline
const std::vector<Atom>&
BulkCrystal::get_basis(void) const
{
  return _basis;
}


inline
const std::vector<Atom>&
BulkCrystal::get_lattice_basis(void) const
{
  return _lattice_basis;
}


inline
const std::vector<Atom>&
BulkCrystal::get_rotated_basis(void) const
{
  return _rotated_basis;
}


inline
const std::string
BulkCrystal::get_lattice_type(void) const
{
  return _lattice_type;
}


inline
const std::vector<double>&
BulkCrystal::get_ortho_lattice_constants(void) const
{
  return _ortho_lattice_constants;
}


inline
const Tensor2Gen& 
BulkCrystal::get_prim_vec(void) const
{
  return _prim_vec;
}

inline
const Tensor2Gen& 
BulkCrystal::get_rotated_prim_vec(void) const
{
  return _rotated_prim_vec;
}


inline
const Tensor2Gen& 
BulkCrystal::get_conv_vec(void) const
{
  return _conv_vec;
}

inline
const Tensor2Gen& 
BulkCrystal::get_rotated_conv_vec(void) const
{
  return _rotated_conv_vec;
}


inline
const Tensor2Gen&
BulkCrystal::get_rotation(void) const
{
  return _rotation;
}

#endif // _BULKCRYSTAL_H_
