// $Id$

#ifndef _TIBEREQSYSTEM_H_
#define _TIBEREQSYSTEM_H_


#include "ModelOptions.h"
#include "IDSet.h"

#include <string>

namespace libMesh
{
class EquationSystems;
class System;
class Elem;
template <typename T> class NumericVector;
template <typename T> class DenseVector;
template <typename T> class DenseMatrix;
}



//! A base class for linear and nonlinear equation systems in TiberCAD
/*!
 * This base class provides some more functionality than the libMesh classes
 * which are useful in TiberCAD.
 */
class TiberEqSystem
{

  public:

    //! The type of system (linear, nonlinear)
    enum SystemType
    {
      UNDEFINED = 0, //!< undefined system type
      LINEAR,        //!< linear system
      NONLINEAR      //!< nonlinear system
    };
    
    //! The type of norms
    enum NormType
    {
      MAX_NORM,  //< the maximum norm
      l2_NORM   //< the l2 norm
    };


    //! Destructor
    virtual ~TiberEqSystem(void) { };


    //! Create a system
    /*!
     * \param[in] es the EquationSystems object
     * \param[in] sysname the name of the new system
     * \param[in] type the type of system (linear, nonlinear)
     * \param[in] options the options for the new system
     * \return a reference to the newly created system
     */

    static TiberEqSystem* create(libMesh::EquationSystems& es,
        const std::string& sysname, SystemType type,
        const ModelOptions& options);


    //! Create a system
    /*!
     * \param[in] es the EquationSystems object
     * \param[in] sysname the name of the new system
     * \param[in] type the type of system (linear, nonlinear)
     * \param[in] options the options for the new system
     * \return a reference to the newly created system
     */
    static TiberEqSystem* create(libMesh::EquationSystems& es,
        const std::string& sysname, const std::string& type,
        const ModelOptions& options);


    //! Set options
    /*!
     * set_options() has to be called after creation of the system.
     *
     * \param options the options as obtained from get_solver_options()
     * in SimulationInterface
     *
     * If the correct options are passed, they will contain a field with
     * key "name" wich is the name of the associated SimulationInterface.
     */
    void set_options(const ModelOptions& options);


    //! Get the options
    const ModelOptions& get_options(void) const;


    //! Get the type of this system
    SystemType get_type(void) const;


    //! Get a pointer to the underlying libmesh System
    /*!
     * In future the return value may be \c NULL if there is no
     * libmesh system type associated
     */
    libMesh::System* get_libmesh_system(void);


    //! Get the solution vector including ghost values
    virtual libMesh::NumericVector<double>& get_solution_vector(void) = 0;


    //! Get the local solution vector without ghost values
    virtual libMesh::NumericVector<double>& get_local_solution_vector(void) = 0;


    //! Set a weight for a given norm
    void set_weight(const libMesh::NumericVector<double>* weight, NormType norm);


    //! Set the set of excluded DoFs
    /*!
     * \param region_ids the region IDs with excluded dofs
     * An empty \c excluded_dofs set resets any previous restrictions
     */
    void set_excluded_dofs(const IDHashSet& excluded_dofs,
        const IDHashSet& interface_dofs = IDHashSet(),
        const std::set<ID>& region_ids = std::set<ID>());

    //! Exclude DoFs from a matrix
    void exclude_dofs(libMesh::DenseMatrix<double>& mat,
        const std::vector<unsigned int>& dof_indices, const libMesh::Elem* elem = NULL);

    //! Exclude DoFs from a matrix
    void exclude_dofs(libMesh::DenseVector<double>& vec,
        const std::vector<unsigned int>& dof_indices, const libMesh::Elem* elem = NULL);


  protected:

    //! Constructor
    TiberEqSystem(void);


    //! Get access to the options
    ModelOptions& get_options(void);


    //! Set the system type
    /*!
     * Call this from derived classes to set the correct system type
     */
    void set_type(SystemType type);


    //! Parse the options
    /*!
     * This method is called from set_options() and can be used to
     * extract some option.
     */
    virtual void parse_options(void) { };


    //! Calculate a norm, considering the weight
    double calculate_norm(libMesh::NumericVector<double>* vec, NormType norm);

    //! Get the weight associated with a given norm
    const libMesh::NumericVector<double>* get_weight(NormType norm) const;



  private:

    //! The options for this system
    ModelOptions _options;


    //! The type of this system (linear, nonlinear)
    SystemType _type;


    //! Weight for the l2 norm
    const libMesh::NumericVector<double>* _l2_weight;


    //! Weight for the l_infty norm
    const libMesh::NumericVector<double>* _linfty_weight;


    //! A set of DoFs which should be excluded from the calculation
    IDHashSet _excluded_dofs;


    //! A set of DoFs which are on the interface between excluded and included regions
    IDHashSet _interface_dofs;


    //! A set of DoFs which should be exc
    std::set<ID> _excluded_region_ids;

};



//
// inline members
//

inline
void
TiberEqSystem::set_options(const ModelOptions& options)
{
  _options += options;
  parse_options();
}


inline
const ModelOptions&
TiberEqSystem::get_options(void) const
{
  return _options;
}


inline
ModelOptions&
TiberEqSystem::get_options(void)
{
  return _options;
}


inline
TiberEqSystem::SystemType
TiberEqSystem::get_type(void) const
{
  return _type;
}


inline
void
TiberEqSystem::set_type(SystemType type)
{
  _type = type;
}


inline
void
TiberEqSystem::set_excluded_dofs(const IDHashSet& exlcuded_dofs,
    const IDHashSet& interface_dofs,
    const std::set<ID>& region_ids)
{
  _excluded_dofs = exlcuded_dofs;
  _interface_dofs = interface_dofs;
  _excluded_region_ids = region_ids;
}


#endif // _TIBEREQSYSTEM_H_
