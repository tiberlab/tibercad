#include "ElementaryPWL.h"

#include "TiberModule.h"

inline
ElementaryPWL::ElementaryPWL(const ModelOptions& options) :
  ElementaryCell(options)
{
}


inline
ElementaryPWL*
ElementaryPWL::create(const ModelOptions& options)
{
  ElementaryPWL* cd = new ElementaryPWL(options);

  return cd;
}



void
ElementaryPWL::do_write_netlist(unsigned int top_node, unsigned int bottom_node,
                                unsigned int& next_free,
                                double area,
                                const libMesh::Point& p,
                                std::ostream& os) const 
{
}
