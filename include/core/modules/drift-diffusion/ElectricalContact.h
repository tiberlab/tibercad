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
class FowlerNordheim;

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


    //! \c true if field emission should be caclulated
    bool has_field_emission(void) const;


    //! Calculate and return the field emission current density
    double calculate_field_emission(double F);


    //! Set the field emission current
    void set_field_emission_current(double J);


    //! Get the field emission current
    double get_field_emission_current(void) const;

    
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


    //! Tells if this is a real current carrying contact
    bool is_real_contact(void) const;


    //! Tells if this is an outer boundary of the simulation domain
    bool is_outer_boundary(void) const;


    //! Tells if this is an internal boundary of the simulation domain
    bool is_internal_boundary(void) const;


    //! Set this contact as outer boundary
    void is_outer_boundary(bool is_outer_boundary);


  protected:

    //! This class is not intended for direct use
    ElectricalContact(void);

    //! Set the boundary condition type for \c variable
    void set_type(DriftDiffusionDefs::Variable variable,
        BCType type);

    //! Get a reference to the Drift-Diffusion properties
    DriftDiffusionProperties& get_material(void);

    //! Get the reference Drift-Diffusion properties
    DriftDiffusionProperties& get_reference_material(void);

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


    //! Set if this is a real current carrying contact
    void is_real_contact(bool is_real_contact);


    //! Determines the reference material for this contact
    void determine_reference_material(void);



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

    //! A pointer to the DriftDiffusionProperties object
    DriftDiffusionProperties *_properties;

    //! The reference material properties
    DriftDiffusionProperties* _reference_prop;


    //! \c true if this is a real (current carrying) contact
    bool _real_contact;


    //! \c true if this is an outer boundary
    bool _is_outer_boundary;


    //! \c true if field emission has to be calculated
    bool _has_field_emission;


    //! The field emission current
    double _field_emission_current;


    //! The field emission model
    FowlerNordheim* _field_emission;

};




//
// inline members
//


inline
ElectricalContact::ElectricalContact(void)
  : _boundary_value(0.0),
    _properties(NULL),
    _real_contact(true),
    _is_outer_boundary(true),
    _has_field_emission(false),
    _field_emission_current(0.0),
    _field_emission(NULL)
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
DriftDiffusionProperties&
ElectricalContact::get_reference_material(void)
{
  return *_reference_prop;
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


inline
void
ElectricalContact::is_real_contact(bool is_real_contact)
{
  _real_contact = is_real_contact;
}


inline
bool
ElectricalContact::is_real_contact(void) const
{
  return _real_contact;
}


inline
bool
ElectricalContact::is_outer_boundary(void) const
{
  return _is_outer_boundary;
}


inline
bool
ElectricalContact::is_internal_boundary(void) const
{
  return !_is_outer_boundary;
}



inline
void
ElectricalContact::is_outer_boundary(bool is_outer_boundary)
{
  _is_outer_boundary = is_outer_boundary;
}



inline
bool
ElectricalContact::has_field_emission(void) const
{
  return _has_field_emission;
}



inline
void
ElectricalContact::set_field_emission_current(double J)
{
  _field_emission_current = J;
}



inline
double
ElectricalContact::get_field_emission_current(void) const
{
  return _field_emission_current;
}


#endif // _ELECTRICALCONTACT_H_
