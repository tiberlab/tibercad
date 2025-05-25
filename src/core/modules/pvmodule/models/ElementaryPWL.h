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

    //! Write the netlist
    virtual void do_write_netlist(unsigned int top_node, unsigned int bottom_node,
                                  unsigned int& next_free,
                                  double area,
                                  const libMesh::Point& p,
                                  std::ostream& os) const override;

  private:

    std::vector<double> _jv_v;

    std::vector<double> _jv_j;

};


#endif // _ELEMENTARYPWL_H_
