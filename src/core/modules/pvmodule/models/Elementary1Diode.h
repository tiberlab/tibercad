#ifndef _ELEMENTARY1DIODE_H_
#define _ELEMENTARY1DIODE_H_

#include "ElementaryCell.h"


class Elementary1Diode : public ElementaryCell
{

  public:

    virtual ~Elementary1Diode(void) override {};

    static Elementary1Diode* create(const ModelOptions& options);


  protected:

    Elementary1Diode(const ModelOptions& options);

    virtual void do_init(void) override;

    //! Write the netlist
    virtual void do_write_netlist(unsigned int top_node, unsigned int bottom_node,
                                  unsigned int& next_free,
                                  double area,
                                  const libMesh::Point& p,
                                  std::ostream& os) const override;

  private:

    //! The series resistance, Ohms*cm^2
    double _rseries = 0.01;

    //! The shunt resistance, Ohms*cm^2
    double _rshunt = 1000;

    //! The photocurrent density, A/cm^2
    double _photocurr = 0.02;

    //! The diode saturation current, A/cm^2
    double _isat = 1e-9;

    //! The diode ideality factor
    double _ideality = 1;
    

};


#endif // _ELEMENTARY1DIODE_H_
