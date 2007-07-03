#ifndef _KSPACE_H_
#define _KSPACE_H_


#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include <vector>
#include "SimulationInterface.h"
#include "tensor.h"
#include <node.h>
class Kspace : public SimulationInterface
{
 public:

  //! Constructor
  Kspace();
  
  //! Destructor
  virtual ~Kspace();


  //!returns reference to kmesh object
  const Mesh& get_k_mesh(void) const; 


  double get_degeneracy_factor(void) {return degeneracy_factor; };

  enum Wedge
  {
    ALL = 0,
    HALF = 1,
    QUARTER =2
  };


 private:

  Kspace(const Kspace& t);
  
  double degeneracy_factor;

  Wedge wedge;

 protected:

  //!integration_order;
  Order integration_order;


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

  //!defines 2D  Brilluoin zone  \f$ k \in [-{\bf K}_1 / 2; {\bf K}_1 / 2) \otimes [-{\bf K}_2 / 2; {\bf K}_2 / 2)   \otimes [-{\bf K}_3 / 2; {\bf K}_3 / 2) \f$
  /*!
    \param k_vector1  Basis k-vector  \f$ {\bf K}_1 \f$ [atom. units]
    \param n - initial number of nodes in direction 1
    \param k_vector2  Basis k-vector  \f$ {\bf K}_2 \f$ [atom. units]
    \param m - initial number of nodes in direction 2
    \param k_vector3  Basis k-vector  \f$ {\bf K}_3 \f$ [atom. units]
    \param k - initial number of nodes in direction 3
  */
  
  void define_k_space(Tensor1 k_vector1, unsigned int n, Tensor1 k_vector2, unsigned int m, Tensor1 k_vector3, unsigned int k);

  
  //!Brilluoin zone
  Mesh* kmesh;


  //!Boundaries of the Brilluoin zone [atomic units]
  double kmin[3], kmax[3];


  //!Dimension of the Brilluoin zone
  short  k_dim;



  //!build k space grid
  void build_k_grid( );


  //!number of nodes in k-domain
  std::vector<unsigned int>  num_nodes;



  //! Rotate mesh
  /*!
    \param mesh  pointer to the mesh
    \RotMatrix transformation matrix (not necessaryly rotation matrix) 
  */
  void rotate_mesh(Mesh* mesh, Tensor2Gen& RotMatrix);


  //! matrix that rotates mesh
  Tensor2Gen transform_matrix;

 
  virtual void 	do_init(void);



  virtual void 	parse_options(void);


};

#endif
