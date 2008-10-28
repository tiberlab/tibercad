// $Id$

#ifndef _DATAOUTPUT_H_
#define _DATAOUTPUT_H_


#include <string>
#include <vector>

class Mesh;


//! A wrapper class for data output
class DataOutput
{

  public:

    //! A type for the known data formats
    enum DataFormat
    {
      UNKNOWN   = 0x0000,       //!< Unknown format
      TECPLOT   = 0x0001,       //!< Tecplot format
      GRACE     = 0x0002,       //!< Xmgrace format
      GNUPLOT   = 0x0004,       //!< GnuPlot format
      VTK       = 0x0008,       //!< Paraview format
      GMSH      = 0x0010,       //!< GMSH format
      GMV       = 0x0020        //!< GMV format
    };


    //! The constructor
    /*!
     * The constructor needs a reference to the mesh and the type of
     * output format.
     *
     * \param mesh the mesh
     * \param format the output file format in string representation
     */
    DataOutput(const Mesh& mesh, const std::string& format);


    //! Write nodal data
    void write_nodal_data(const std::string& filename,
        const std::vector<double>& data,
        const std::vector<std::string>& legend);


    //! Write cell data
    void write_cell_data(const std::string& filename,
        const std::vector<double>& data,
        const std::vector<std::string>& legend);


    //! Get data format for a given data format name
    static DataFormat tell_data_format(const std::string& format);


    //! Set the output directory
    void set_output_directory(const std::string& output_dir);



  private:

    //! The mesh
    const Mesh* _mesh;

    //! The file format
    unsigned int _format;

    //! The output directory
    std::string _output_dir;

};


#endif // _DATAOUTPUT_H_
