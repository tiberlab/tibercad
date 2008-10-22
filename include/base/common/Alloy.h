// $Id$

#ifndef _ALLOY_H_
#define _ALLOY_H_


#include "Material.h"


//! Description of  an  alloy  material.
/*!
 * An alloy is thought to be of the form \f$A_xB_{1-x}C\f$.
 * 
 * This class contains the properties of the constituent Materials A and B.
 */
class Alloy : public Material
{

  public:
    
    //! Destructor
    /*!
     * Deletes all \c PhysicalProperties objects
     */
    virtual ~Alloy(void);
   

    //! Create a material with name \c name
    static Alloy* create(const std::string& name);


    //! Return the name of component material A
    const std::string& get_name_A(void) const;
    
    //! Return the name of component material B
    const std::string& get_name_B(void) const;


  protected:

    //! Construct an  alloy material 
    Alloy(const std::string& name);


    //! \copydoc Material::do_init()
    virtual void do_init(void);


  private:

    //! Molar fraction x of the \c Alloy
    /*!
     *  This is the \f$x\f$ in \f$A_xB_{1-x}C\f$
     */
    double _molar_fraction;


    //! The component A
    Material* _mat_A;

    
    //! The component B
    Material* _mat_B;

};



//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------



inline
const std::string&
Alloy::get_name_A(void) const
{
  return _mat_A->get_name();
}


inline
const std::string&
Alloy::get_name_B(void) const
{
  return _mat_B->get_name();
}


#endif /* _ALLOY_H_ */
