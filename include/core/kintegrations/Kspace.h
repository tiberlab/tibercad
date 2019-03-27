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


  enum Wedge
  {
    ALL = 0,
    HALF = 1,
    QUARTER = 2,
    EIGHTH = 3
  };

  //! Rotate mesh
  void rotate_mesh(void);

  void inv_rotate_mesh(void);

  //calculates angular between 2D k space basis (y-z) vectors to different between quadratic and hexagonal k space
  void define_type_of_k_space();

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
     HEXAGONAL      //< 2/3D hexagonal

   };

   //!build k space grid
   void build_k_grid( );

   //!defines 1D  Brilluoin zone \f$ k \in [-{\bf K}/2; {\bf K}/2) \f$
   /*!
    \param k_vector Basis k-vector  \f$ \bf K \f$ [atom. units]
    \param n - initial number of nodes
    */
   void define_k_space(Tensor1 k_vector, unsigned int n);

   //!defines 2D  Brilluoin zone  \f$ k \in [-{\bf K}_1 / 2; {\bf K}_1 / 2) \otimes [-{\bf K}_2 / 2; {\bf K}_2 / 2) \f$
   /*!
    \param k_vector1  Basis k-vector  \f$ {\bf K}_1 \f$ [atom. units]
    \param n - initial number of nodes in direction 1
    \param k_vector2  Basis k-vector  \f$ {\bf K}_2 \f$ [atom. units]
    \param m - initial number of nodes in direction 2
    */
   void define_k_space(Tensor1 k_vector1,unsigned int n,  Tensor1 k_vector2, unsigned int m);

   //!defines 3D  Brilluoin zone  \f$ k \in [-{\bf K}_1 / 2; {\bf K}_1 / 2) \otimes [-{\bf K}_2 / 2; {\bf K}_2 / 2)   \otimes [-{\bf K}_3 / 2; {\bf K}_3 / 2) \f$
   /*!
    \param k_vector1  Basis k-vector  \f$ {\bf K}_1 \f$ [atom. units]
    \param n - initial number of nodes in direction 1
    \param k_vector2  Basis k-vector  \f$ {\bf K}_2 \f$ [atom. units]
    \param m - initial number of nodes in direction 2
    \param k_vector3  Basis k-vector  \f$ {\bf K}_3 \f$ [atom. units]
    \param k - initial number of nodes in direction 3
    */
   void define_k_space(Tensor1 k_vector1, unsigned int n, Tensor1 k_vector2, 
                       unsigned int m, Tensor1 k_vector3, unsigned int k);

   void define_k_path(void);
    
   void tokenize(const std::string& str,
                std::vector<std::string>& tokens,
                const std::string& delimiters = "-");

   ModelOptions mod_opt;

   double degeneracy_factor;

   Wedge wedge;

   //! mesh_order
   libMeshEnums::Order _mesh_order;


   //! Brilluoin zone
   libMesh::Mesh* kmesh;

   //! Boundaries of the Brilluoin zone
   double kmin[3], kmax[3];

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

   //! Generic identifier vector
   /*!
    * This vector is used to identify the different basis
    * vector according to the k-space symmetry.
    */
   std::vector<unsigned int> identification;

   //! The symmetry class of the k-space
   Symmetry k_space_symmetry;

   //! True if this is a path in k-space
   bool k_path;

   //! The communicator for the k-space
   libMesh::Parallel::Communicator kspace_comm;

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
