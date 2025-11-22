// $Id$


#ifndef _XMONITOR_H_
#define _XMONITOR_H_


#include "petsc.h"

#include <string>

//! A general X monitor
class XMonitor
{

  public:

    //! Create an X monitor
    static XMonitor* create(const std::string& title);

    
    //! Destructor
    ~XMonitor(void);


    //! Set axis labels
    void set_axis_labels(const std::string& x, const std::string& y);
    
    //! Draw a point
    void draw_point(double x, double y);
    

  private:

    //! Constructor
    explicit XMonitor(const std::string& title);


    //! PetscDraw object
    PetscDraw _draw;

    //! PetscDrawLG object
    PetscDrawLG _lg;

    //! View ports
    PetscDrawViewPorts* _ports;

    //! The axis handle
    PetscDrawAxis _axis;

    //! The title
    std::string _title;

};

//
// inline methods
//



inline
XMonitor*
XMonitor::create(const std::string& title)
{
  XMonitor* mon = NULL;
  try
  {
    mon = new XMonitor(title);
  }
  catch (...)
  {
  }

  return mon;
}


inline
void
XMonitor::set_axis_labels(const std::string& x, const std::string& y)
{
  if (_axis != NULL)
    PetscDrawAxisSetLabels(_axis, _title.c_str(), x.c_str(), y.c_str());
}


inline
void
XMonitor::draw_point(double x, double y)
{
  PetscDrawLGAddPoint(_lg, &x, &y);
  PetscDrawLGDraw(_lg);
}


#endif // _XMONITOR_H_
