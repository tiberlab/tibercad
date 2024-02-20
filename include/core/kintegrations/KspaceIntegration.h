#ifndef _KSPACEINTEGRATION_H_
#define _KSPACEINTEGRATION_H_


#include "TiberModelObject.h"
#include "Kspace.h"
#include "TypeDefs.h"
#include "libMeshDefs.h"

#include "libmesh/error_vector.h"
#include "libmesh/elem.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/dof_object.h"
#include "libmesh/parallel.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>


#include <fstream>
#include <iomanip>
#include <map>



//! Definition of DofField that maps a dof index to a double.
typedef std::vector<double> DofField;
// TO DO Change in NumericVector:
//typedef NumericVector<double> DofField;

/*! 
 * General class used to perform Kspace Integrations. 
 * It is used with Kspace and KspaceIntegrationTemplate
 * The class can perform integrations over a generic Kspace
 * of any field. For generality the field is a general vector<double>.
 *
 * Adaptive Mesh refinements over the Kspace are possible using libmesh 
 * for this calculate_for_k_point( ) must return into "double estimator"
 * a meaningful error estimator scalar at the given k_point.   
 */
class KspaceIntegration : public TiberModelObject 
{
 public:

  typedef Elem KElem;
  typedef std::map <const KElem*, double> KMeshToIntegratedValue;



  //! Creator for KspaceIntegration objects
  /*!
   * \param hook the pointer to the object which will calculate quantities
   *        at each k-point
   * \param callback the member function of the above class that does the calculation
   * \param opts options to control the integration
   */
  template <class T>
  static KspaceIntegration* create(T* hook,
      void (T::*callback)(const Point&, DofField&, double&),
      const ModelOptions& opts,
      const libMesh::Parallel::Communicator& k_comm);

  template <class T>
  static KspaceIntegration* create(T* hook,
      void (T::*callback)(const Point&, const Point&, DofField&, double&),
      const ModelOptions& opts,
      const libMesh::Parallel::Communicator& k_comm);


  virtual ~KspaceIntegration();
  
  void init(void);
  
  void solve(void);

  const DofField& get_solution(void) const;
 
  void get_solution(DofField& density) const;

  const Kspace* get_k_space(void) const;

  Kspace* get_k_space(void);

  unsigned int get_k_space_dimension(void) const;

  //!creates a parallel communicator by splitting the Device communicator (larger than mesh_comm)
  //!all nodes with same id (color) of the mesh_communicator have to compute the same k-point
  static void create_communicator(const libMesh::Parallel::Communicator& device_comm,
                                  const libMesh::Parallel::Communicator& mesh_comm,
                                  libMesh::Parallel::Communicator& communicator);

  //!returns the k points and weights
  std::map<libMesh::Point, double>
  get_kpoints_and_weights(void);
  
  //!returns the kpoints and weights of the kspace as a shared global quantity among the processors.
  std::map<DofField, double>
  broadcast_kpoints(libMesh::Parallel::Communicator& comm, std::map<Point, double> k_points);

  //!distributes the k-grid across the processors and for each process returns the local indices
  std::vector<int>
  distribute_kpoints(const libMesh::Parallel::Communicator& comm, std::map<DofField, double> k_points);

  //! builds the map local k points - weights
  std::map<DofField, double>
  get_local_kpoints(std::map<DofField, double> global_kpoints, std::vector<int> local_k_indices);


 protected:

  KspaceIntegration(const ModelOptions& options, 
                    const libMesh::Parallel::Communicator& k_comm);

  virtual void do_init(void);

  virtual void do_solve(void);

   //!put data into opt
  virtual void parse_options(void);
 
  //!calculates the quantity performing mesh refinement of k-space, if required.
  virtual void calculate_convergent_density(void);

  //!calculates density that is necessary for eack k-point and a number that will be used for refinement 
  virtual void calculate_for_k_point(const Point& k_point,
      const Point& refpoint,
      DofField& density,
      double& estimator) = 0;

  
  virtual double estimate_error(void);

  //!calculate density for a particular k-grid
  virtual void calculate_density();

  //!calculates integrated quantity distribution over k space
  //virtual std::vector<double>  get_density_in_k_space(void)  const;

  
 private:


  //!options for charge density
  struct options
  {
    //! if true user provides Brillouin zone size. Otherwise, it is calculated by the program
    bool k_domain_user_input;

    //! if true, program will refine the Brilluoin zone adaptively
    bool k_domain_refinement;

    //! if true, all the cells in the k-space are refined
    bool uniform_refinement;

    //! fraction of the elements to be refined
    double refine_fraction;

    //! stop refinement if \f$ ||\rho_{i+1} - \rho_i||/||rho_i|| < \epsilon \f$
    double relative_accuracy;

    //! maximum level for k space refinement
    unsigned int maximum_ref_level;

    //! degeneracy factor
    unsigned int degeneracy;

    //!perflog output
    bool log_output;

    //! normalization volume
    double normalization_volume;
  };

  //! KIntegration options
  options opt; 


  //!result after integration: old map
  DofField old_density; 

  //!integration_order;
  libMeshEnums::QuadratureType quadrature_type;
   

  //!integration_order;
  libMeshEnums::Order integration_order;
  

  //!fem_order;
  libMeshEnums::Order fem_order;
  

  //!result after integration: maps the real space mesh to the integrated quantity
  DofField real_space_density;

  //! temporary field arrays
  DofField dens_at_k_elem;

  DofField dens_at_k_point;

  //! kspace must be built in the derived classes (??)
  Kspace* _kspace;

  //!maps k-space to a real space density (which is a map between real space elements and density)   
  std::map< const KElem*, DofField > density_at_k;


  //!map from Elem in the k-grid to an integrated quantity used for the refinement criterion
  KMeshToIntegratedValue error_estimator;

  //!is used for k-space output
  std::string additional_name_suffix;

  //!how many elementsin k-space has to be done
  unsigned int count_elements() const;
   
  //!calculate  volumes of the elements
  void calculate_volumes(void) const;
   
  //!estimates error for mesh refinement KellyErrorEstimator is called
  void estimate_error_for_refinement(libMesh::ErrorVector& error);

  libMesh::Parallel::Communicator kspace_comm;

};


#include "KspaceIntegrationTemplate.h"

template <class T>
KspaceIntegration*
KspaceIntegration::create(T* hook,
    void (T::*callback)(const Point&, DofField&, double&),
    const ModelOptions& opts,
    const libMesh::Parallel::Communicator& k_comm)
{
  KspaceIntegrationTemplate<T>* templ = new KspaceIntegrationTemplate<T>(hook, opts, k_comm);
  templ->_callback = callback;
  return templ;
}

template <class T>
KspaceIntegration*
KspaceIntegration::create(T* hook,
    void (T::*callback)(const Point&, const Point&, DofField&, double&),
    const ModelOptions& opts,
    const libMesh::Parallel::Communicator& k_comm)
{
  KspaceIntegrationTemplate<T>* templ = new KspaceIntegrationTemplate<T>(hook, opts, k_comm);
  templ->_callback2 = callback;
  return templ;
}



inline
void
KspaceIntegration::init(void)
{
  do_init();
}

inline
void
KspaceIntegration::solve(void)
{
  do_solve();
}

inline
const Kspace*
KspaceIntegration::get_k_space(void) const
{
  return(_kspace);
}

inline
Kspace*
KspaceIntegration::get_k_space(void)
{
  return(_kspace);
}


inline
const DofField&
KspaceIntegration::get_solution(void) const 
{
  return real_space_density;
}

inline
void
KspaceIntegration::get_solution(DofField& density) const
{
  density = real_space_density;
}

#endif
