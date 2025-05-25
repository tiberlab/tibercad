
#ifndef _ELEMENTARYCELL_H_
#define _ELEMENTARYCELL_H_

#include "PhysicalModel.h"

#include <ostream>

class Point;

// Base class for elementary cell models
class ElementaryCell : public PhysicalModel
{

  public:

    virtual ~ElementaryCell(void) {};

    //! Write the netlist
    /*!
     * This method writes the elementary cell netlist into the
     * circuit file.
     *
     * \param top_node the index of the top layer node
     * \param bottom_node the index of the bottom layer node
     * \param next_free the next free circuit node, will be updated during
     *        the call
     * \param area the cell area in cm^2
     * \param p the coordinates of the elementary cell
     * \param os the stream to write to
     */
    void write_netlist(unsigned int top_node, unsigned int bottom_node,
                       unsigned int& next_free,
                       double area,
                       const libMesh::Point& p,
                       std::ostream& os) const
    { do_write_netlist(top_node, bottom_node, next_free, area, p, os); };


  protected:

    ElementaryCell(const ModelOptions& options);


    //! Write the netlist
    /*!
     * This method writes the elementary cell netlist into the
     * circuit file.
     *
     * \param top_node the index of the top layer node
     * \param bottom_node the index of the bottom layer node
     * \param next_free the next free circuit node, will be updated during
     *        the call
     * \param area the cell area in cm^2
     * \param p the coordinates of the elementary cell
     * \param os the stream to write to
     */

    virtual void do_write_netlist(unsigned int top_node, unsigned int bottom_node,
                                  unsigned int& next_free,
                                  double area,
                                  const libMesh::Point& p,
                                  std::ostream& os) const = 0;

  private:


};


inline
ElementaryCell::ElementaryCell(const ModelOptions& options) :
  PhysicalModel(options)
{
}



#endif // _ELEMENTARYCELL_H_
