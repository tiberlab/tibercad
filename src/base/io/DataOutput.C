// $Id$

#include "DataOutput.h"

#include "GMVIO_cell.h"
#include "tecplot_IO_cell.h"
#include "gnuplot_io.h"
#include "GraceIO.h"
#include "VTKIO.h"


#include "mesh.h"


DataOutput::DataOutput(const Mesh& mesh, const std::string& format)
  : _mesh(&mesh),
    _format(format)
{
}




void
DataOutput::write_nodal_data(const std::string& filename,
    const std::vector<Number>& data,
    const std::vector<std::string>& legend)
{
  // default is GMV
  
  if (_format == "ise")
    TecplotIO(*_mesh).write_nodal_data(filename + ".plt", data, legend);
  else if (_format == "grace")
    GraceIO(*_mesh).write_nodal_data(filename + ".dat", data, legend);
  else if (_format == "gnuplot")
    GnuPlotIO(*_mesh).write_nodal_data(filename + ".dat", data, legend);
  else if (_format == "vtk")
    VTKIO(*_mesh).write_nodal_data(filename + ".vtk", data, legend);
  else
    GMVIO(*_mesh).write_nodal_data(filename + ".gmv", data, legend);
}





void
DataOutput::write_cell_data(const std::string& filename,
    const std::vector<Number>& data,
    const std::vector<std::string>& legend)
{
  // default is GMV
  
  if (_format == "ise")
    TecplotIO_cell(*_mesh).write_cell_data(filename + ".plt", data, legend);
  else if (_format == "grace")
    GraceIO(*_mesh).write_elemental_data(filename + ".dat", data, legend);
  else if (_format == "gnuplot")
    std::cout << "GnuPlot does not currently support cell data." << std::endl;
  else if (_format == "vtk")
    VTKIO(*_mesh).write_elemental_data(filename + ".vtk", data, legend);
  else
    GMVIO_cell(*_mesh).write_ascii_cell_data(filename + ".gmv", data, legend);
}


