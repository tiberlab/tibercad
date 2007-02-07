// $Id$

#ifndef _SWEEPABLE_H_
#define _SWEEPABLE_H_

//! Interface for classes which contain a sweepable variable
class Sweepable
{

  public:

    //! Empty constructor
    Sweepable(void);

    //! Destructor
    virtual ~Sweepable(void);
    
    //! Set a new value
    void set_new_value(double value);

    //! Get the current value
    double get_current_value(void) const;

    //! Get the old value
    double get_old_value(void) const;


  private:

    //! The current value of the sweep variable
    double _value;

    //! The old value of the sweep variable
    double _old_value;

};



//
// inline members
// 

inline 
Sweepable::Sweepable(void)
  : _value(0.0),
    _old_value(0.0)
{
}


inline
Sweepable::~Sweepable(void)
{
}


inline
void
Sweepable::set_new_value(double value)
{
  _old_value = _value;
  _value = value;
}


inline
double
Sweepable::get_current_value(void) const
{
  return _value;
}


inline
double
Sweepable::get_old_value(void) const
{
  return _old_value;
}



#endif // _SWEEPABLE_H_
