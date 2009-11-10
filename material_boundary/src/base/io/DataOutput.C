// $Id$

#include "DataOutput.h"
#include "Utils.h"

#include "GMVIO_cell.h"
#include "tecplot_IO_cell.h"
#include "gnuplot_io.h"
#include "GraceIO.h"
#include "VTKIO.h"


#include "mesh.h"
#include "gmsh_io.h"


DataOutput::DataOutput(const MeshBase& mesh, const std::string& format)
  : _mesh(&mesh),
    _format(0x0000)
{
  std::vector<std::string> formats;
  Utils::extract_vector(format, formats);

  for (int i = 0; i < formats.size(); i++)
    _format |= tell_data_format(formats[i]);

}


DataOutput::DataFormat
DataOutput::tell_data_format(const std::string& format)
{
  DataFormat df = UNKNOWN;

  if ((format == "ise") || (format == "tecplot"))
    df = TECPLOT;
  else if (format == "grace")
    df = GRACE;
  else if ((format == "vtk") || (format == "paraview"))
    df = VTK;
  else if (format == "gmv")
    df = GMV;
  else if (format == "gmsh")
    df = GMSH;
  else if (format == "gnuplot")
    df = GNUPLOT;


  return df;
}



void
DataOutput::write_nodal_data(const std::string& filename,
    const std::vector<Number>& data,
    const std::vector<std::string>& legend)
{
  std::string file(_output_dir);
  file += "/" + filename;

  if (_format & TECPLOT)
    TecplotIO(*_mesh).write_nodal_data(file + ".plt", data, legend);
  if (_format & GRACE)
    GraceIO(*_mesh).write_nodal_data(file + ".dat", data, legend);
  if (_format & GNUPLOT)
    GnuPlotIO(*_mesh).write_nodal_data(file + ".dat", data, legend);
  if (_format & VTK)
    TiberVTKIO(*_mesh).write_nodal_data(file + ".vtk", data, legend);
  if (_format & GMSH)
    GmshIO(*_mesh).write_nodal_data(file + ".msh", data, legend);
  if (_format & GMV)
    GMVIO(*_mesh).write_nodal_data(file + ".gmv", data, legend);
}





void
DataOutput::write_cell_data(const std::string& filename,
    const std::vector<Number>& data,
    const std::vector<std::string>& legend)
{
  std::string file(_output_dir);
  file += "/" + filename;

  if (_format & TECPLOT)
    TecplotIO_cell(*_mesh).write_cell_data(file + ".plt", data, legend);
  if (_format & GRACE)
    GraceIO(*_mesh).write_elemental_data(file + ".dat", data, legend);
  if (_format & GNUPLOT)
    std::cout << "GnuPlot does not currently support cell data." << std::endl;
  if (_format & VTK)
    TiberVTKIO(*_mesh).write_elemental_data(file + ".vtk", data, legend);
  if (_format & GMV)
    GMVIO_cell(*_mesh).write_ascii_cell_data(file + ".gmv", data, legend);
}




void
DataOutput::set_output_directory(const std::string& output_dir)
{
  _output_dir = output_dir;
}



void
DataOutput::set_mesh(const MeshBase& mesh)
{
  _mesh = &mesh;
}
