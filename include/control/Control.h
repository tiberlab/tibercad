// $Id$

#ifndef _CONTROL_H_
#define _CONTROL_H_

#include "ModelErrorException.h"
#include "InitFailedException.h"
#include "DatabaseException.h"
#include "SolveFailedException.h"
#include "TypeDefs.h"
#include "IDSet.h"

#include <map>
#include <list>
#include <vector>
#include <set>
#include <iostream>
#include <cassert>

class Device;
class Database;
class SimulationInterface;
class SimulationEnvironment;

//! The control module of TiberCAD
/*!
 * This class is responsible for the setup of the device and
 * the flow control of the simulation
 */
class Control
{
  // Some typedefs that are needed in the public section
  private:

    //! Type for the list of simulations
    typedef std::map<const std::string, SimulationInterface*> SimulationMap;


    //! Type for the list of simulation environments
    typedef std::map<SimulationInterface*,
            SimulationEnvironment*> EnvironmentMap;


  public:

    //! An iterator to iterate over all simulations
    class simulation_iterator
    {
      public:
        simulation_iterator(void) {};
        simulation_iterator(const simulation_iterator& it) : _it(it._it) {};
        simulation_iterator(const SimulationMap::iterator& it) : _it(it) {};

        simulation_iterator& operator=(const simulation_iterator& other)
        {
          _it = other._it;
          return *this;
        };

        simulation_iterator& operator++(void) { ++_it; return *this; };

        simulation_iterator& operator--(void) { --_it; return *this; };

        SimulationInterface* operator*(void) { return _it->second; };

        bool operator==(const simulation_iterator& other)
        {
          return (_it == other._it);
        };

        bool operator!=(const simulation_iterator& other)
        {
          return (_it != other._it);
        };

      private:
        SimulationMap::iterator _it;
    };


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


    //! Find a simulation with name \c name
    /*!
     * \param name the name to look for
     * \return a pointer to the simulation if found, \c NULL otherwise
     *
     * \c name can be one of the following:
     * \li the user defined name of a simulation
     * \li the identifier of the simulation as used for creation
     * \li the empty string
     *
     * In the second case, the first simulation of this type will be
     * returned. In the third case, the first of all simulations will
     * be returned.
     *
     */
    SimulationInterface* find_simulation(const std::string& name) const;


    //! Get the directory where to put output files
    const std::string& get_output_dir(void) const;


    //! Get the complete suffix for the output filenames
    const std::string& get_filename_suffix(void) const;


    //! Clear the suffix for the output filenames
    void clear_filename_suffix(void);


    //! Append something to the suffix for the output filenames
    /*!
     * The filename suffix will be appended to all output files which
     * contain plot data.
     * The suffix itself will be prepended by a '_'
     */
    void append_to_filename_suffix(const std::string& suffix);


    //! Prepend something to the suffix for the output filenames
    /*!
     * The filename suffix will be appended to all output files which
     * contain plot data.
     * The suffix itself will be prepended by a '_'
     */
    void prepend_to_filename_suffix(const std::string& suffix);


    //! Delete the first output filename suffix part
    void drop_first_filename_suffix(void);


    //! Delete the last output filename suffix part
    void drop_last_filename_suffix(void);



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


    //! Mark all environments as unprepared
    void invalidate_environments(void);


    //! Get the iterator to the first simulation
    simulation_iterator simulations_begin(void);


    //! Get the past-the-end iterator for the simulations
    simulation_iterator simulations_end(void);


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


    //! A list of all simulations we control
    SimulationMap _simulations;


    //! A list of all simulation environments we control
    /*!
     * \note Not every simulation necessarily has an associated
     * environment!
     */
    EnvironmentMap _simulation_environments;


    //! The list of simulations to be solved
    std::vector<std::string> _solve_list;


    //! The directory where to put output
    std::string _outputdir;


    //! The list from which the filename suffix gets constructed
    std::list<std::string> _filename_suffix;


    //! The filename suffix string
    /*!
     * This is only a helper variable
     */
    mutable std::string _filename_suffix_str;


    //! The output format
    /*!
     * see get_output_format() for a detailed description
     */
    std::string _output_format;

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
void
Control::clear_filename_suffix(void)
{
  _filename_suffix.clear();
}



inline
void
Control::append_to_filename_suffix(const std::string& suffix)
{
  if (suffix.size() != 0)
    _filename_suffix.push_back(suffix);
}


inline
void
Control::prepend_to_filename_suffix(const std::string& suffix)
{
  if (suffix.size() != 0)
    _filename_suffix.push_front(suffix);
}


inline
void
Control::drop_last_filename_suffix(void)
{
  if (_filename_suffix.size() != 0)
    _filename_suffix.pop_back();
}


inline
void
Control::drop_first_filename_suffix(void)
{
  if (_filename_suffix.size() != 0)
    _filename_suffix.pop_front();
}



inline
const std::string&
Control::get_output_format(void) const
{
  return _output_format;
}



inline
Control::simulation_iterator
Control::simulations_begin(void)
{
  return simulation_iterator(_simulations.begin());
}



inline
Control::simulation_iterator
Control::simulations_end(void)
{
  return simulation_iterator(_simulations.end());
}



#endif // _CONTROL_H_
