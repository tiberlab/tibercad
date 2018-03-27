// $Id$

#ifndef _TIBERCAD_H_
#define _TIBERCAD_H_

#include "tiber_dll.h"
//#include "TiberMath.h"

#include <string>
#include <list>

#include "parallel.h"

class Control;

namespace libMesh
{
class LibMeshInit;
}


//! Entry point and generic stuff for TiberCAD
class TiberCad
{

  public:

    //! The statistics to be used
    enum Statistics
    {
      BOLTZMANN,     /*!< Boltzmann statistics */
      FERMIDIRAC,    /*!< Fermi-Dirac statistics */
      BOSEEINSTEIN   /*!< Bose-Einstein statistics */
    };


    //! Types of symmetries
    enum Symmetry
    {
        NONE,         /*!< no special symmetry */
        CYLINDRICAL   /*!< cylinder symmetry */
    };


    //! Constructor
    /*
     * \param mpi_comm the MPI communicator we should use
     */
    TiberCad(MPI_Comm mpi_comm);


    //! Destructor
    ~TiberCad(void);



    //! An initialization routine
    /*!
     * This routine calls init() of libmesh and other libraries, if needed.
     */
    void init(const std::string& inputfile);


    //! Start simulations
    void run(void);



    //! The full TiberCAD version string
    static std::string version_string(bool include_svn_release = false);


    //! The architecture as string
    static std::string arch_string(void);


    //! The major version number
    static int major_version(void);


    //! The minor version number
    static int minor_version(void);


    //! The subminor version number
    static int subminor_version(void);


    //! The software revision version number
    static int software_revision(void);


    //! Get the directory where to put output files
    static const std::string& get_output_dir(void);


    //! Get the complete suffix for the output filenames
    static std::string get_filename_suffix(void);


    //! Clear the suffix for the output filenames
    static void clear_filename_suffix(void);


    //! Append something to the suffix for the output filenames
    /*!
     * The filename suffix will be appended to all output files which
     * contain plot data.
     * The suffix itself will be prepended by a '_'
     */
    static void append_to_filename_suffix(const std::string& suffix);


    //! Prepend something to the suffix for the output filenames
    /*!
     * The filename suffix will be appended to all output files which
     * contain plot data.
     * The suffix itself will be prepended by a '_'
     */
    static void prepend_to_filename_suffix(const std::string& suffix);


    //! Delete the first output filename suffix part
    static void drop_first_filename_suffix(void);


    //! Delete the last output filename suffix part
    static void drop_last_filename_suffix(void);


    //! Get the output format
    /*!
     * \return a string that identifies the type of output files
     * to generate
     *
     * Currently the following formats are supported:
     * \li \c gmv for GMV
     * \li \c ise for Tecplot
     * \li \c gnu for GnuPlot
     */
    static const std::string& get_output_format(void);


    //! Get our MPI Communicator
    static libMesh::Parallel::Communicator& get_mpi_comm(void);



  private:

    //! A cleanup routine
    /*!
     * This routine calls close() of libmesh and other libraries, if needed
     */
    void cleanup(void);


    //! The TiberCAD major version
    static const int _MajorVersion;


    //! The TiberCAD minor version
    static const int _MinorVersion;


    //! The TiberCAD subminor version
    static const int _SubMinorVersion;


    //! The TiberCAD subversion release number
    static const int _SvnRevision;


    //! A counter to assure that there is only one instance of this class
    static unsigned int _object_counter;


    //! The installation root directory
    /*!
     * This is read from the environment
     */
    static std::string _tiberroot;

    //! The flow control object
    static Control* _control;

    //! The libmesh entry point
    libMesh::LibMeshInit* _libmeshinit;


    //! The list from which the filename suffix gets constructed
    static std::list<std::string> _filename_suffix;


    //! This is our MPI communicator
    static libMesh::Parallel::Communicator _mpi_comm;


};


#endif // _TIBERCAD_H_
