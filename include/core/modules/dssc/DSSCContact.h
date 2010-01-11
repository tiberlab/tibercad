// $Id$

#ifndef _DSSCCONTACT_H_
#define _DSSCCONTACT_H_

#include "BoundaryProperties.h"


class DSSCContact : public BoundaryProperties
{

  public:

    //! Constructor
    DSSCContact(const ModelOptions& options);

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

    bool is_open_circuit(void) const;


  protected:

    /*! \copydoc BoundaryProperties::do_init() */
    virtual void do_init(void);



  private:

    //! The boundary value (eg. applied voltage)
    double _boundary_value;

    bool _cathode;

    double _Ioc;
    double _Idark;

    double _I3oc;
    double _I3dark;

    double _j0;

    double _beta;

    double _current;

    static bool _open_circuit;


    //! Set the open circuit potential and densities
    void set_OC_values(double Ioc, double Idark, double I3oc, double I3dark);

    void calculate_current(double I, double I3);



};


//
// inline members
//

inline
DSSCContact::DSSCContact(const ModelOptions& options)
  // open circuit value
  : BoundaryProperties(options),
    _boundary_value(1e10),
    _cathode(false),
    _j0(0.1),
    _beta(0.78),
    _current(0.0)
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
  //std::cerr << "Ioc = "  <<Ioc << std::endl
  //  << "I3oc = "  <<I3oc << std::endl;
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
    calculate_current(I, I3);
}




inline
double
DSSCContact::get_current(void) const
{
  return _current;
}



inline
bool
DSSCContact::is_open_circuit(void) const
{
  return _open_circuit;
}


#endif // _DSSCCONTACT_H_
