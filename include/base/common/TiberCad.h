// This header is part of the public tibercad API

#ifndef _TIBERCAD_H_
#define _TIBERCAD_H_

#include "tiber_dll.h"

#include <string>
#include <list>

#include "libmesh/parallel.h"

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

    //! A cleanup routine
    /*!
     * This routine calls close() of libmesh and other libraries, if needed
     */
    void cleanup(void);


    //! Start simulations
    void run(void);


    //! Get the current global time (i.e. simulation time coordinate)
    static double get_global_time(void);


    //! The full TiberCAD version string
    static std::string version_string(bool include_compilation_date = false);


    //! The architecture as string
    static std::string arch_string(void);


    //! The compilation date
    static const std::string& compilation_date(void);


    //! The compilation system
    static const std::string& compilation_system(void);


    //! Last modification date
    static const std::string& last_modification(void);


    //! The software revision version number
    static const std::string& software_revision(void);


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


    //! Get our MPI Communicator
    static libMesh::Parallel::Communicator& get_mpi_comm(void);



  private:

    //! The TiberCAD subversion release number
    static const std::string _git_revision;


    //! Last modified date
    static const std::string _git_modified;


    //! Compilation date
    static const std::string _compilation_date;


    //! Compilation system
    static const std::string _compilation_system;


    //! A counter to assure that there is only one instance of this class
    static unsigned int _object_counter;

    //! Control variable to check whether we are initialized
    static bool _is_initialized;


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
