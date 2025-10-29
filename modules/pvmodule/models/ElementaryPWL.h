#ifndef _ELEMENTARYPWL_H_
#define _ELEMENTARYPWL_H_

#include "ElementaryCell.h"


class ElementaryPWL : public ElementaryCell
{

  public:

    virtual ~ElementaryPWL(void) override {};

    static ElementaryPWL* create(const ModelOptions& options);


  protected:

    ElementaryPWL(const ModelOptions& options);

    virtual void do_init(void) override;

    //! Write the netlist
    virtual void do_write_netlist(unsigned int top_node, unsigned int bottom_node,
                                  unsigned int& next_free,
                                  double area,
                                  const libMesh::Elem* elem,
                                  const libMesh::Point& p,
                                  std::ostream& os) const override;

  private:

    //! The voltage data points
    std::vector<double> _jv_v;

    //! The current density data points
    std::vector<double> _jv_j;

};


#endif // _ELEMENTARYPWL_H_
