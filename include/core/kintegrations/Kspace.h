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

class Kspace
{
 public:

  //! Constructor
  Kspace( const ModelOptions& mod_opt);

  //! copy constr
  Kspace( const Kspace& kspace);

  //! Destructor
  virtual ~Kspace();


  //!returns reference to kmesh object
  //Mesh* get_k_mesh(void);

  const Mesh* get_k_mesh(void);

  short mesh_dimension(void);

  double get_degeneracy_factor(void);

  //libMeshEnums::Order get_integration_order(void);
  //libMeshEnums::QuadratureType get_quadrature_type();

  enum Wedge
  {
    ALL = 0,
    HALF = 1,
    QUARTER =2,
    EIGHTH = 3
  };

  //! Rotate mesh
  /*!
    \param mesh  pointer to the mesh
    \RotMatrix transformation matrix (not necessaryly rotation matrix)
  */
  void rotate_mesh(void);

  void inv_rotate_mesh(void);

 protected:

   virtual void  do_init(void) throw (InitFailedException);

   virtual void  parse_options(void);


 private:

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

   //!mesh_order;
   libMeshEnums::Order mesh_order;


   //!Brilluoin zone
   Mesh* kmesh;

   //!Boundaries of the Brilluoin zone [atomic units]
   double kmin[3], kmax[3];

   //!Dimension of the k_space
   unsigned int  k_dim;

   //!Dimension of the Brilluoin zone
   unsigned int  k_space_dim;

   //!number of nodes in k-domain
   std::vector<unsigned int>  num_nodes;

   //! matrix that rotates mesh
   Tensor2Gen transform_matrix;

   bool k_path;


};

//---------------------------------------------------------------------//
inline  double Kspace::get_degeneracy_factor(void)
{
    return  degeneracy_factor;
}

//inline libMeshEnums::Order Kspace::get_integration_order(void)
//{
//   return integration_order;
//}

//inline libMeshEnums::QuadratureType Kspace::get_quadrature_type(void)
//{
//   return quadrature_type;
//}
//inline Mesh* Kspace::get_k_mesh()
//{
 // return(kmesh);
//}

inline const Mesh* Kspace::get_k_mesh()
{
  return kmesh;
}

inline short Kspace::mesh_dimension(void)
{
  return k_dim;
}

#endif
