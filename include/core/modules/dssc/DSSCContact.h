// $Id$

#ifndef _DSSCCONTACT_H_
#define _DSSCCONTACT_H_

#include "BoundaryProperties.h"
#include "Variable.h"


class DSSCContact : public BoundaryProperties, public Variable
{

  public:

    //! Constructor
    DSSCContact(void);

    //! Destructor
    ~DSSCContact(void) { };

    //! Create a contact
    static DSSCContact* create(const std::string& name,
        const ModelOptions& options = ModelOptions());


    //! Get the contact potential
    double get_potential(void) const;


    //! Set the potential
    void set_potential(double potential);


    bool& is_cathode(void);


  protected:

    /*! \copydoc BoundaryProperties::do_init() */
    virtual void do_init(void);


    /*! \copydoc Variable::set_variable_value() */
    virtual void set_variable_value(double value, ID id = 0);


    /*! \copydoc Variable::get_variable_value() */
    virtual double get_variable_value(ID id = 0);



  private:
    
    //! The boundary value (eg. applied voltage)
    double _boundary_value;
  
    bool _cathode;
    
};


//
// inline members
// 

inline
DSSCContact::DSSCContact(void)
  : _boundary_value(0.0),
    _cathode(false)
{
}


inline
void
DSSCContact::set_potential(double potential)
{
  _boundary_value = potential;
}


inline
double
DSSCContact::get_potential(void) const
{
  return _boundary_value;
}



inline
void
DSSCContact::set_variable_value(double value, ID id)
{
  ignore_unused_variable(id);
  set_potential(value);
}


inline
double
DSSCContact::get_variable_value(ID id)
{
  ignore_unused_variable(id);
  return get_potential();
}



inline
bool&
DSSCContact::is_cathode(void)
{
  return _cathode;
}

#endif // _DSSCCONTACT_H_
