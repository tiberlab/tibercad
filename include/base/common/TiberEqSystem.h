// $Id$

#ifndef _TIBEREQSYSTEM_H_
#define _TIBEREQSYSTEM_H_


#include "ModelOptions.h"
#include "IDSet.h"

#include <string>

class EquationSystems;
template <typename T> class NumericVector;
template <typename T> class DenseVector;
template <typename T> class DenseMatrix;


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
      l2_NORM,   //< the l2 norm
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
    static TiberEqSystem* create(EquationSystems& es,
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
    static TiberEqSystem* create(EquationSystems& es,
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


    //! Get the solution vector
    virtual NumericVector<double>& get_solution_vector(void) = 0;


    //! Set a weight for a given norm
    void set_weight(const NumericVector<double>* weight, NormType norm);


    //! Set the set of excluded DoFs
    void set_excluded_dofs(const IDHashSet& exlcuded_dofs);

    //! Exclude DoFs from a matrix
    void exclude_dofs(DenseMatrix<double>& mat, const std::vector<unsigned int>& dof_indices);

    //! Exclude DoFs from a matrix
    void exclude_dofs(DenseVector<double>& vec, const std::vector<unsigned int>& dof_indices);


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
    double calculate_norm(NumericVector<double>* vec, NormType norm);

    //! Get the weight associated with a given norm
    const NumericVector<double>* get_weight(NormType norm) const;



  private:

    //! The options for this system
    ModelOptions _options;


    //! The type of this system (linear, nonlinear)
    SystemType _type;


    //! Weight for the l2 norm
    const NumericVector<double>* _l2_weight;


    //! Weight for the l_infty norm
    const NumericVector<double>* _linfty_weight;


    //! A set of DoFs which should be excluded from the calculation
    IDHashSet _excluded_dofs;

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
TiberEqSystem::set_excluded_dofs(const IDHashSet& exlcuded_dofs)
{
  _excluded_dofs = exlcuded_dofs;
}


#endif // _TIBEREQSYSTEM_H_
