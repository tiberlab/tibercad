#ifndef _KSPACE_H_
#define _KSPACE_H_

#include "enum_order.h"
#include "enum_quadrature_type.h"
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include <vector>
#include "tensor.h"
#include <node.h>
#include "ModelOptions.h"
#include "InitFailedException.h"
#include "parallel.h"

class Kspace
{
 public:

  //! Constructor
  Kspace( const ModelOptions& mod_opt, const libMesh::Parallel::Communicator& comm_in);

  //! copy constr
  Kspace( const Kspace& kspace);

  //! Destructor
  virtual ~Kspace();


  //! returns reference to kmesh object
  libMesh::MeshBase* get_k_mesh(void);

  //! returns reference to kmesh object
  const libMesh::MeshBase* get_k_mesh(void) const;

  //! Get the dimensionality of the k-space
  unsigned int dimension(void) const;

  //! Get the degeneracy
  double get_degeneracy_factor(void);

  //! Tells whether this is a k-path or not
  bool is_k_path(void) const;

  //! Transform a k-point to relative coordinates
  void inverse_transform(libMesh::Point& p) const;

  //! Get coordinates for a high symmetry point
  /*!
   * The coordinates will be given in the reciprocal basis as
   * calculated for the k-space.
   * The provided string should be the standard notation for
   * symmetry points for the given crystal type. In case an invalid
   * point is provided, an exception will be thrown.
   */
  libMesh::Point get_symmetry_point(const std::string& name) const;

  //! Transform a point according to the rotation matrix
  void transform_point(libMesh::Point& p) const;

  //! Rotate mesh
  void rotate_mesh(void);

  void inv_rotate_mesh(void);

 protected:

   virtual void  do_init(void);

   virtual void  parse_options(void);


 private:

   //! Known symmetry types of the Brillouin zone
   enum Symmetry
   {
     GAMMA     = 0, //< Gamma point only
     LINEAR,        //< 1D is always linear
     QUADRATIC,     //< 2D quadratic
     RECTANGULAR,   //< 2D rectangular
     HEXAGONAL,     //< 2/3D hexagonal
     FCC,
     BCC,
     CUBIC,
     TETRAGONAL,
     ORTHORHOMBIC

   };

   //! Detect the symmetry class of the kspace
   void find_k_space_symmetry(void);

   //!build k space grid
   void build_k_grid(void);

   //!defines 1D  Brilluoin zone \f$ k \in [-{\bf K}/2; {\bf K}/2) \f$
   /*!
    \param k_vector Basis k-vector  \f$ \bf K \f$ [atom. units]
    */
   void define_k_space(Tensor1 k_vector);

   //!defines 2D  Brilluoin zone  \f$ k \in [-{\bf K}_1 / 2; {\bf K}_1 / 2) \otimes [-{\bf K}_2 / 2; {\bf K}_2 / 2) \f$
   /*!
    \param k_vector1  Basis k-vector  \f$ {\bf K}_1 \f$ [atom. units]
    \param k_vector2  Basis k-vector  \f$ {\bf K}_2 \f$ [atom. units]
    */
   void define_k_space(Tensor1 k_vector1, Tensor1 k_vector2);

   //!defines 3D  Brilluoin zone  \f$ k \in [-{\bf K}_1 / 2; {\bf K}_1 / 2) \otimes [-{\bf K}_2 / 2; {\bf K}_2 / 2)   \otimes [-{\bf K}_3 / 2; {\bf K}_3 / 2) \f$
   /*!
    \param k_vector1  Basis k-vector  \f$ {\bf K}_1 \f$ [atom. units]
    \param k_vector2  Basis k-vector  \f$ {\bf K}_2 \f$ [atom. units]
    \param k_vector3  Basis k-vector  \f$ {\bf K}_3 \f$ [atom. units]
    */
   void define_k_space(Tensor1 k_vector1, Tensor1 k_vector2,
                       Tensor1 k_vector3);

   //! Define a linear path in k space for dispersions
   void define_k_path(void);
    
   //! The options
   ModelOptions mod_opt;

   //! Degeneracy factor for reduced wedge of the BZ
   double degeneracy_factor;

   //! mesh_order
   libMeshEnums::Order _mesh_order;


   //! Brilluoin zone
   libMesh::Mesh* kmesh;

   //! Dimension of the k_space
   unsigned int  k_space_dim;

   //! number of nodes in k-domain
   std::vector<unsigned int>  num_nodes;

   //! matrix that rotates mesh
   Tensor2Gen transform_matrix;

   //! k-space basis vector 1
   libMesh::Point k_basis_vector1;

   //! k-space basis vector 2
   libMesh::Point k_basis_vector2;

   //! k-space basis vector 3
   libMesh::Point k_basis_vector3;

   //! Identify k basis vector b1
   /*!
    * This and \c b2 and \c b3 are used to identify the basis
    * vectors for the different crystal classes, since the calculated
    * basis depends on the orientation of the system. This allows to easily
    * calculate symmetry points using standard coordinates afterwards.
    * C indexing is used (0, 1, 2)
    */
   unsigned int b1;

   //! Identify k basis vector b1
   unsigned int b2;

   //! Identify k basis vector b1
   unsigned int b3;

   //! The maximum k to use (relative to k space basis)
   libMesh::RealVectorValue k_max;

   //! The symmetry class of the k-space
   Symmetry k_space_symmetry;

   //! True if this is a path in k-space
   bool k_path;

   //! The communicator for the k-space
   libMesh::Parallel::Communicator kspace_comm;

   //! The symmetries as strings
   static std::vector<std::string> symmetry_names;

};

//---------------------------------------------------------------------//
inline  double Kspace::get_degeneracy_factor(void)
{
    return  degeneracy_factor;
}


inline bool Kspace::is_k_path(void) const
{
  return(k_path);
}

inline libMesh::MeshBase* Kspace::get_k_mesh()
{
  return(kmesh);
}

inline const libMesh::MeshBase* Kspace::get_k_mesh() const
{
  return(kmesh);
}

inline unsigned int Kspace::dimension(void) const
{
  return(k_space_dim);
}

#endif
