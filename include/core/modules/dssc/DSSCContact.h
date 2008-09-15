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

    void set_values(double V, double I, double I3);


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

    double _Voc;
    double _Ioc;
    double _I3oc;

    double _current;

    bool _open_circuit;
    

    //! Set the open circuit potential and densities
    void set_OC_values(double Voc, double Ioc, double I3oc);

    void calculate_current(void);

};


//
// inline members
// 

inline
DSSCContact::DSSCContact(void)
  : _boundary_value(0.0),
    _cathode(false),
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
DSSCContact::set_OC_values(double Voc, double Ioc, double I3oc)
{
  _Voc = Voc;
  _Ioc = Ioc;
  _I3oc = I3oc;
  std::cerr << "Voc = " << Voc << std::endl
    << "Ioc = "  <<Ioc << std::endl
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
DSSCContact::set_values(double V, double I, double I3)
{
  if (_open_circuit)
  {
    set_OC_values(V, I, I3);
    _current = 0.0;
  }
  else
    calculate_current();
}




inline
double
DSSCContact::get_current(void) const
{
  return _current;
}



#endif // _DSSCCONTACT_H_
