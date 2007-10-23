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


  private:

    //! The mesh
    const Mesh* _mesh;

    //! The file format
    std::string _format;

};


#endif // _DATAOUTPUT_H_
