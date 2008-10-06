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

    //! Get the contact current
    double get_current(void) const;


    //! Set the potential
    void set_potential(double potential);


    bool& is_cathode(void);

    void set_values(double I, double Idark, double I3, double I3dark);


    void set_open_circuit(bool open_circuit = true);



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

    double _Ioc;
    double _Idark;

    double _I3oc;
    double _I3dark;

    double _current;
   
    bool _open_circuit;
    

    //! Set the open circuit potential and densities
    void set_OC_values(double Ioc, double Idark, double I3oc, double I3dark);

    void calculate_current(void);


 
};


//
// inline members
// 

inline
DSSCContact::DSSCContact(void)
  : _boundary_value(0.0),
    _cathode(false),
    _current(0.0),
    _open_circuit(true)
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
  set_open_circuit(false);
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

inline
void
DSSCContact::set_OC_values(double Ioc, double Idark, double I3oc, double I3dark)
{
  _Ioc = Ioc;
  _Idark = Idark;
  _I3oc = I3oc;
  _I3dark = I3dark;
  std::cerr << "Ioc = "  <<Ioc << std::endl
    << "I3oc = "  <<I3oc << std::endl;
}


inline
void
DSSCContact::set_open_circuit(bool open_circuit)
{
  _open_circuit = open_circuit;
}


inline
void
DSSCContact::set_values(double I, double Idark, double I3, double I3dark)
{
  if (_open_circuit)
  {
    set_OC_values(I, Idark, I3, I3dark);
    _current = 0.0;
  }
  else
    calculate_current();
}




inline
double
DSSCContact::get_current(void) const
{
  return _boundary_value;
  //return _current;
}



#endif // _DSSCCONTACT_H_
