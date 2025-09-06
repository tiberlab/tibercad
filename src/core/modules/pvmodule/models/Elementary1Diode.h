#ifndef _ELEMENTARY1DIODE_H_
#define _ELEMENTARY1DIODE_H_

#include "ElementaryCell.h"

class Photocurrent;
class DegradationModel;

class Elementary1Diode : public ElementaryCell
{

  public:

    virtual ~Elementary1Diode(void) override {};

    static Elementary1Diode* create(const ModelOptions& options);


  protected:

    Elementary1Diode(const ModelOptions& options);

    virtual void do_init(void) override;
    
    virtual void prepare_submodels(void) override;

    //! Write the netlist to the provided stream
    virtual void do_write_netlist(unsigned int top_node, unsigned int bottom_node,
                                  unsigned int& next_free,
                                  double area,
                                  const libMesh::Elem* elem,
                                  const libMesh::Point& p,
                                  std::ostream& os) const override;

  private:

    //! The series resistance, Ohms*cm^2
    double _rseries = 1;

    //! The shunt resistance, Ohms*cm^2
    double _rshunt = 10000;

    //! The photocurrent density, A/cm^2
    double _photocurr = 0.02;

    //! The diode saturation current, A/cm^2
    double _isat = 1e-11;

    //! The diode ideality factor
    double _ideality = 2;
    
    //! We can have a submodel for the photocurrent
    Photocurrent* _photocurr_model = nullptr;
    
    //! We can have a submodel for the cell degradation
    /*!
     * In future we might allow for several model, acting
     * on different parameters based on different quantities
     */
    DegradationModel* _degradation_model = nullptr;

};


#endif // _ELEMENTARY1DIODE_H_
