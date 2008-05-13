// $Id$

#ifndef _TIBEREQSYSTEM_H_
#define _TIBEREQSYSTEM_H_


#include "ModelOptions.h"

#include <string>

class EquationSystems;

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
     * \param[in] sysname the name of teh new system
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
    void set_options(const ModelOptions& options);


    //! Get the options
    const ModelOptions& get_options(void) const;


    //! Get the type of this system
    SystemType get_type(void) const;



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


  private:

    //! The options for this system
    ModelOptions _options;


    //! The type of this system (linear, nonlinear)
    SystemType _type;

};



//
// inline members
//

inline
void
TiberEqSystem::set_options(const ModelOptions& options)
{
  _options += options;
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


#endif // _TIBEREQSYSTEM_H_
