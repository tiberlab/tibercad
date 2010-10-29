// $Id$

#ifndef _CONTROL_H_
#define _CONTROL_H_

#include "ModelErrorException.h"
#include "InitFailedException.h"
#include "DatabaseException.h"
#include "SolveFailedException.h"
#include "TypeDefs.h"
#include "IDSet.h"
#include "tiber_dll.h"

#include <map>
#include <list>
#include <vector>
#include <set>
#include <iostream>
#include <cassert>

class Device;
class Database;

//! The control module of TiberCAD
/*!
 * This class is responsible for the setup of the device and
 * the flow control of the simulation
 */
class TBDLLOCAL Control
{

  public:


    //! The constructor
    /*!
     * \param inputfile the input file to be used
     */
    Control(void);


    //! The destructor
    ~Control(void);


    //! Set the input file
    void set_inputfile(const std::string& inputfile);


    //! Initialize the device and the simulations
    /*!
     * This method calls create_device(), create_materials() and
     * setup_models()
     */
    void init(void) throw (InitFailedException,
        ModelErrorException, DatabaseException);


    //! Runs the simulation
    /*!
     * Calls solve() of the simulations according to the rules given
     * in the input file
     */
    void run_simulation(void) throw (SolveFailedException);


    //! Plots the results of all simulation models
    void plot_all(void);


    //! Get a reference to the device
    Device& get_device(void);



    //! Get the directory where to put output files
    const std::string& get_output_dir(void) const;



    //! Get the output format
    /*!
     * \return a string that identifies the type ouf output files
     * to generate
     *
     * Currently the following formats are supported:
     * \li \c gmv for GMV
     * \li \c ise for Tecplot
     * \li \c gnu for GnuPlot
     */
    const std::string& get_output_format(void) const;



  private:

    //! A signal handler
    class SignalHandler
    {

      public:

        typedef struct sigaction SigAction;

        static void set_control(Control* ctrl);

        static void activate_sigint(void);
        static void deactivate_sigint(void);

      private:

        SignalHandler(void);

        static Control* _ctrl;

        static SigAction _old_int_action;

        static void sigint(int sig);
    };


    //! The inputfile
    std::string _inputfile;


    //! The device we control
    Device* _device;


    //! The database we use
    Database* _database;


    //! The list of simulations to be solved
    std::vector<std::string> _solve_list;


    //! The directory where to put output
    std::string _outputdir;


    //! The output format
    /*!
     * see get_output_format() for a detailed description
     */
    std::string _output_format;
  
    //! Rotate from xy to xz plane in case of cylindrical symmetry.
    bool _xy2xz_rotation;

    //! Setup global simulation options
    void setup_globals(void);

    //! Create the device
    /*!
     * The creation of a device involves the following function calls
     * in the order as given here:
     * \code
     * create_device();
     * create_materials();
     * create_atomistic_structures();
     * setup_clusters();
     * setup_models();
     * \endcode
     */
    void create_device(void);


    //! Create all the materials
    void create_materials(void);


    //! Create all atomistic structures
    void create_atomistic_structures(void);


    //! Create and setup the models
    void setup_models(void) throw (InitFailedException, ModelErrorException);


    //! Setup region clusters
    /*!
     * A region cluster contains different physical regions of the mesh
     * which possibly overlap with the material regions described in the
     * \c Region sections
     */
    void setup_clusters(void);


    //! Extract physical regions from a string
    void extract_physical_regions(const std::string& str, IDSet& ids);
};


//
// inline members
//


inline
void
Control::set_inputfile(const std::string& inputfile)
{
  _inputfile = inputfile;
}


inline
Device&
Control::get_device(void)
{
  return *_device;
}



inline
const std::string&
Control::get_output_dir(void) const
{
  return _outputdir;
}




inline
const std::string&
Control::get_output_format(void) const
{
  return _output_format;
}





#endif // _CONTROL_H_
