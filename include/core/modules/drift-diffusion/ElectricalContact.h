// $Id$

#ifndef _ELECTRICALCONTACT_H_
#define _ELECTRICALCONTACT_H_


#include "DriftDiffusionDefs.h"
#include "BoundaryProperties.h"
#include "ModelOptions.h"
#include "Variable.h"

// C++ includes
#include <string>

class DriftDiffusionProperties;

//! Base class for any kind of electrical contact models
/*!
 * Derive from this class to define contact models. Examples can be found
 * in OhmicContact and SchottkyContact.
 *
 * This class is derived also from Variable to be able to make
 * a voltage sweep.
 */
class ElectricalContact : public BoundaryProperties, public Variable
{
  public:

    //! Empty destructor
    virtual ~ElectricalContact(void) {};

    //! Create an electrical contact model
    static ElectricalContact* create(const std::string& name,
        const ModelOptions& options = ModelOptions());

    
    //! The type of boundary condition
    enum BCType
    {
      DIRICHLET, /*!< Dirichlet type BC */
      NEUMANN,   /*!< von Neumann type BC */
      MIXED,     /*!< mixed type BC */
      //! pinning of electro-chemical potential with respect to 
      //! the conduction band edge
      PINNING    
    };

    //! Set the simulation voltage for this contact
    void set_simulation_voltage(double voltage);

    //! Get the simulation voltage for this contact
    double get_simulation_voltage(void) const;

    //! Get the inner voltage for this contact
    /*!
     * This takes surface resistance into account
     */
    double get_inner_voltage(void) const;
    
    //! Set the material
    /*!
     * The material reference is needed to compute parameters of the boundary
     * condition
     */
    void set_material(DriftDiffusionProperties *properties);


    //! Set the particle fluxes normal to the surface
    void set_normal_fluxes(double jn, double jp);

    
    //! Get the type of boundary condition for \c variable
    BCType get_type(DriftDiffusionDefs::Variable variable) const;

    //! Get the normal derivative for \c variable
    /*!
     * Starting from general boundary conditions
     * \f[\alpha u + \beta\frac{\partial u}{\partial\vec{n}}=\gamma\f]
     * we write the normal derivative as
     * \f[\frac{\partial u}{\partial\vec{n}}=c - au\f]
     * where \f$c=\frac{\gamma}{\beta}\f$ and \f$a=\frac{\alpha}{\beta}\f$
     *
     * NOTE: For the electric potential we assume \c c to be a surface charge
     * density in units of cm^-2
     *
     * \param[in] variable the variable for which coefficients should be returned
     * \param[out] a the coefficient a
     * \param[out] c the coefficient c
     */
    virtual void get_normal_derivative(DriftDiffusionDefs::Variable variable,
        double& a, double& c);


    //! Get the derivatives of the parameters for the normal derivative
    virtual void get_derivatives_of_normal_derivative(
        DriftDiffusionDefs::Variable variable,
        std::vector<double>& da, std::vector<double>& dc);
        

    //! Get the boundary value for \c variable
    virtual double get_boundary_value(DriftDiffusionDefs::Variable variable);

    //! Set a BC to homogeneous von Neumann type
    /*!
     * This is only useful if on has a Dirichlet type boundary condition
     * and wants to change it into homogeneous von Neumann.
     * It can be used to have different BCs for electron and hole
     * electro-chemical potentials in small devices, or to change to zero
     * electric field BC
     */
    void set_zero_derivative_bc(DriftDiffusionDefs::Variable variable);


  protected:

    //! This class is not intended for direct use
    ElectricalContact(void);

    //! Set the boundary condition type for \c variable
    void set_type(DriftDiffusionDefs::Variable variable,
        BCType type);

    //! Get a reference to the Drift-Diffusion properties
    DriftDiffusionProperties& get_material(void);

    /*! \copydoc BoundaryProperties::do_init() */
    virtual void do_init(void);


    /*! \copydoc Variable::set_variable_value() */
    virtual void set_variable_value(double value, ID id = 0);


    /*! \copydoc Variable::get_variable_value() */
    virtual double get_variable_value(ID id = 0);


    //! Get the normal electron flux
    double get_normal_electron_flux(void) const;


    //! Get the normal hole flux
    double get_normal_hole_flux(void) const;


    //! Get the contact voltage drop
    /*!
     * A positive value means: the inner voltage is lower than
     * the outer contact voltage
     */
    double get_contact_voltage_drop() const;


  private:

    BCType _potential_type;
    BCType _fermie_type;
    BCType _fermih_type;


    //! The boundary value (eg. applied voltage)
    double _boundary_value;
    

    //! The normal electron flux
    double _jn;
    

    //! The normal hole flux
    double _jp;


    //! The contact surface resistance
    double _surfres;

    // A pointer to the DriftDiffusionProperties object
    DriftDiffusionProperties *_properties;
};




//
// inline members
//


inline
ElectricalContact::ElectricalContact(void)
  : _boundary_value(0.0),
    _properties(NULL)
{
}

inline
void
ElectricalContact::set_simulation_voltage(double voltage)
{
  _boundary_value = voltage;
}


inline
double
ElectricalContact::get_simulation_voltage(void) const
{
  return _boundary_value;
}


inline
double
ElectricalContact::get_inner_voltage(void) const
{
  double v = get_simulation_voltage();
  if (_surfres > 1e-12)
    v -= get_contact_voltage_drop();
  
  return v;
}



inline
void
ElectricalContact::set_material(DriftDiffusionProperties *properties)
{
  _properties = properties;
}


inline
void
ElectricalContact::set_normal_fluxes(double jn, double jp)
{
  _jn = jn;
  _jp = jp;
}


inline
double
ElectricalContact::get_normal_electron_flux(void) const
{
  return _jn;
}


inline
double
ElectricalContact::get_normal_hole_flux(void) const
{
  return _jp;
}



inline
DriftDiffusionProperties&
ElectricalContact::get_material(void)
{
  return *_properties;
}

inline
ElectricalContact::BCType
ElectricalContact::get_type(DriftDiffusionDefs::Variable variable) const
{
  BCType type;

  switch (variable)
  {
    case DriftDiffusionDefs::POTENTIAL:
      type = _potential_type;
      break;
    case DriftDiffusionDefs::FERMIE:
      type = _fermie_type;
      break;
    case DriftDiffusionDefs::FERMIH:
      type = _fermih_type;
      break;
  }
  
  return type;
}


inline
void
ElectricalContact::set_type(DriftDiffusionDefs::Variable variable, BCType type)
{
  switch (variable)
  {
    case DriftDiffusionDefs::POTENTIAL:
      _potential_type = type;
      break;
    case DriftDiffusionDefs::FERMIE:
      _fermie_type = type;
      break;
    case DriftDiffusionDefs::FERMIH:
      _fermih_type = type;
      break;
  }
}

inline
void
ElectricalContact::get_normal_derivative(DriftDiffusionDefs::Variable variable,
        double& a, double& c)
{
  ignore_unused_variable(variable);
  a = 0.0;
  c = 0.0;
}

inline
double
ElectricalContact::get_boundary_value(DriftDiffusionDefs::Variable variable)
{
  ignore_unused_variable(variable);
  return 0.0;
}

inline
void
ElectricalContact::set_zero_derivative_bc(DriftDiffusionDefs::Variable variable)
{
  set_type(variable, ElectricalContact::NEUMANN);
}



inline
void
ElectricalContact::set_variable_value(double value, ID id)
{
  ignore_unused_variable(id);
  set_simulation_voltage(value);
}


inline
double
ElectricalContact::get_variable_value(ID id)
{
  ignore_unused_variable(id);
  return get_simulation_voltage();
}


inline
void
ElectricalContact::get_derivatives_of_normal_derivative(
    DriftDiffusionDefs::Variable variable,
    std::vector<double>& da, std::vector<double>& dc)
{
  da = std::vector<double>(3, 0.0);
  dc = std::vector<double>(3, 0.0);
}



#endif // _ELECTRICALCONTACT_H_
