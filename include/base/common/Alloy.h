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


    //! Build virtual cristal approximation of a parameter including bowing
    /*!
     * In a ternary or quaternary compound semiconductor
     * \f$Q = A_xB_{1-x}C\f$ the value of a
     * material parameter can (in the virtual crystal approximation) be
     * calculated as
     * \f[\alpha_Q = x\alpha_{AC} + (1-x)\alpha_{BC} - bx(1-x)\f]
     * where \em b is called bowing parameter and describes deviation
     * from the nonlinear behaviour.
     *
     * \param ac the value for material \f$AC\f$
     * \param bc the value for material \f$BC\f$
     * \param x the molar fraction of \f$AC\f$
     * \param bowing the bowing parameter
     */
    //static double calculate_VCA_parameter(double ac, double bc,
    //    double x, double bowing = 0.0);


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
    std::string _name_A;

    //! The component B
    std::string _name_B;

    //! The models for component A
    ModelMap _models_A;

    //! The models for component B
    ModelMap _models_B;

    //! The rotated crystal of A
    RotatedCrystal* _cryst_A;

    //! The rotated crystal of A
    RotatedCrystal* _cryst_B;
    
};



//--------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------


//inline
//double
//Alloy::calculate_VCA_parameter(double ac, double bc,
//    double x, double bowing)
//{
//  return bc + (ac - bc) * x - bowing * x * (1 - x);
//}


inline
const std::string&
Alloy::get_name_A(void) const
{
  return _name_A;
}


inline
const std::string&
Alloy::get_name_B(void) const
{
  return _name_B;
}


#endif /* _ALLOY_H_ */
