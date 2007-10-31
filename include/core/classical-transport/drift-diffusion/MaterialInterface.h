// $Id$

#ifndef _MATERIALINTERFACE_H_
#define _MATERIALINTERFACE_H_


#include "ElectricalContact.h"


//! A simple class to impose a current density on a boundary
/*!
 * For now one can only impose electron current density, and in the
 * input file it is called voltage
 */
class MaterialInterface : public ElectricalContact
{

  public:

    //! The constructor
    MaterialInterface(void);
    
    //! The destructor
    ~MaterialInterface(void) {};

    //! Create an object
    static MaterialInterface* create(void);

    //! Get normal derivative
    virtual void get_normal_derivative(DriftDiffusionDefs::Variable variable,
        double& a, double& c);

    //! Get derivatives of the normal derivative
    virtual void get_derivatives_of_normal_derivative(
        DriftDiffusionDefs::Variable variable,
        std::vector<double>& da, std::vector<double>& dc);


  protected:

    /*! \copydoc ElectricalContact::do_init() */
    virtual void do_init(void);


  private:

    //! The surface state sheet density
    double _Ns;

    //! The energy of the states with respect to the conduction band edge
    double _Es;

    //! Multiplicity
    double _g_factor;

};


//
// inline members
// 

inline
MaterialInterface*
MaterialInterface::create(void)
{
  return new MaterialInterface();
}



inline
MaterialInterface::MaterialInterface(void)
  : _Ns(0.0),
    _Es(0.2),
    _g_factor(2)
{
  set_type(DriftDiffusionDefs::POTENTIAL, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIE, ElectricalContact::NEUMANN);
  set_type(DriftDiffusionDefs::FERMIH, ElectricalContact::NEUMANN);
}



#endif // _MATERIALINTERFACE_H_
