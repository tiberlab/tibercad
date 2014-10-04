// $Id$

#ifndef _DATAOUTPUT_H_
#define _DATAOUTPUT_H_

#include "SolutionDescriptor.h"

#include <string>
#include <vector>
#include <map>

class MeshBase;


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
      //GNUPLOT   = 0x0004,       //!< GnuPlot format
      VTK       = 0x0008,       //!< Paraview format
      //GMSH      = 0x0010,       //!< GMSH format
      GMV       = 0x0020        //!< GMV format
    };


    //! Constructor
    DataOutput(void);

    //! Destructor
    virtual ~DataOutput(void);


    //! The constructor
    /*!
     * The constructor needs a reference to the mesh and the type of
     * output format.
     *
     * \param mesh the mesh
     * \param format the output file format in string representation
     */
    DataOutput(const MeshBase& mesh, const std::string& format);


    //! Create a writer for a given format
    static DataOutput* create(const std::string& format);


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

    //! Set the filename
    void set_filename(const std::string& filename);

    //! Set to binary output
    void set_binary(void);

    //! Set to ASCII output
    void set_ascii(void);


    //! Set the mesh
    void set_mesh(const MeshBase& mesh);


    //! Set the data to print
    /*!
     * \param data the map containing all data of a certain zone
     * \param zone the zone ID (e.g. mesh subdomain)
     *
     * Nodal data in \c data has to be in the same order as when doing
     * a loop over the elements and then over nodes of the element.
     *
     * Only a reference to the \c data object is kept internally, i.e. the
     * provided data container must persist until after calling \c write().
     */
    void set_data(const std::map<SolutionDescriptor,
        std::vector<double> >& data, ID zone);


    //! Write out everything
    /*!
     * When \c force \c = \c true, file is written even if there
     * is no data.
     */
    void write(bool force = false);



  protected:

    //! A typedef for the data container
    typedef std::map<SolutionDescriptor, std::vector<double> > DataMap;

    //! Get a reference to the mesh
    const MeshBase& get_mesh(void) const;

    //! Check if this writer has a mesh assigned
    bool has_mesh(void) const;

    //! Get the output directory
    const std::string& get_output_directory(void) const;

    //! Get the filename
    /*!
     * \return the full path to the output file, excluding the suffix
     */
    std::string get_filename(void) const;

    //! Check if it is a binary file
    bool is_binary(void) const;

    //! Check if it is an ascii file
    bool is_ascii(void) const;

    //! Check if a zone ID has data
    /*!
     * If no argument is given, the method checks if there is
     * any zone data available.
     */
    bool has_data(ID zone = INVALID_ID);

    //! Get the zone data
    const DataMap& get_zone_data(ID zone);

    //! The actual writer
    virtual void do_write(bool force) {}; // = 0;



  private:

    //! The mesh
    const MeshBase* _mesh;

    //! The file format
    unsigned int _format;

    //! The output directory
    std::string _output_dir;

    //! The filename
    std::string _filename;

    //! \c true if it is binary file
    bool _is_binary;

    //! A map collecting pointers to all zone data
    std::map<ID, const DataMap*> _data;


};


//
// inline members
//

inline
bool
DataOutput::has_mesh(void) const
{
  return (_mesh != NULL ? true : false);
}

inline
bool
DataOutput::is_binary(void) const
{
  return _is_binary;
}

inline
bool
DataOutput::is_ascii(void) const
{
  return !is_binary();
}


#endif // _DATAOUTPUT_H_
