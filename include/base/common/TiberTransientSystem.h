// $Id: TiberTransientSystem.h 4433 2017-07-04 09:37:54Z maufder $

#ifndef _TIBERTRANSIENTSYSTEM_H_
#define _TIBERTRANSIENTSYSTEM_H_

#include "TiberEqSystem.h"
#include "libMeshDefs.h"

#include "libmesh/transient_system.h"




/*!
 * \brief First order transient system in tibercad
 *
 */
template <class Base>
class TiberTransientSystem : public TiberEqSystem, public libMesh::TransientSystem<Base>
{

  public:

    //! Constructor
    TiberTransientSystem(libMesh::EquationSystems& es,
        const std::string& name,
        const unsigned int number);


    //! Destructor
    virtual ~TiberTransientSystem(void);


    //! Create a linear system
    /*!
     * \param es the EquationSystems object where the new system will be added
     * \param sysname the name of the new system
     * \param options the options for this system
     * \return a pointer to the newly created system
     */
    static TiberTransientSystem* create(libMesh::EquationSystems& es,
        const std::string& sysname, const ModelOptions& options);


    /*! \copydoc ImplicitSystem::solve() */
    virtual void solve(void) override;


    /*! \copydoc ImplicitSystem::system_type() */
    virtual std::string system_type(void) const override;


    /*! \copydoc System:user_initialization() */
    virtual void user_initialization(void) override;


    //! Get the solution vector including ghost values
    virtual libMesh::NumericVector<double>& get_solution_vector(void) override;


    //! Get the local solution vector without ghost values
    virtual libMesh::NumericVector<double>& get_local_solution_vector(void) override;


    //! Set the target time
    void set_target_time(double time);


  private:

    typedef libMesh::TransientSystem<Base> parent_type;

    //! The next target time
    double _target_time = 0.;

};
    

typedef TiberTransientSystem<libMesh::LinearImplicitSystem> TiberTransientLinSystem;

#endif // _TIBERTRANSIENTSYSTEM_H_
