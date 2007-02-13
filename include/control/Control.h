// $Id$

#ifndef _CONTROL_H_
#define _CONTROL_H_

#include "ModelErrorException.h"
#include "InitFailedException.h"
#include "SolveFailedException.h"

#include <map>
#include <set>

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

  public:

    //! The constructor
    /*!
     * \param inputfile the input file to be used
     */
    Control(const std::string& inputfile);

    
    //! The destructor
    ~Control(void);

    
    //! Initialize the device and the simulations
    /*!
     * This method calls create_device(), create_materials() and
     * setup_models()
     */
    void init(void) throw (InitFailedException);

    
    //! Runs the simulation
    /*!
     * Calls solve() of the simulations according to the rules given
     * in the input file
     */
    void run_simulation(void) throw (SolveFailedException);
    

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


    //! Get the variables to plot
    const std::set<std::string>& get_plotvariables(void) const;


    //! Get the directory where to put output files
    const std::string& get_output_dir(void) const;


    //! Get the suffix for the output filenames
    const std::string& get_filename_suffix(void) const;


    //! Clear the suffix for the output filenames
    void clear_filename_suffix(void);


    //! Set the suffix for the output filenames
    /*!
     * The filename suffix will be appended to all output files which
     * contain plot data.
     * The suffix itself will be prepended by a '_'
     */
    void set_filename_suffix(const std::string& suffix);


    //! Append something to the suffix for the output filenames
    /*!
     * The filename suffix will be appended to all output files which
     * contain plot data.
     * The suffix itself will be prepended by a '_'
     */
    void append_to_filename_suffix(const std::string& suffix);



  private:

    //! Type for the list of simulations
    typedef std::map<const std::string, SimulationInterface*> SimulationMap;

    
    //! Type for the list of simulation environments
    typedef std::map<SimulationInterface*,
            SimulationEnvironment*> EnvironmentMap;

    
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
     * \note {Not every simulation necessarily has an associated
     * environment! }
     */
    EnvironmentMap _simulation_environments;


    //! The directory where to put output
    std::string _outputdir;


    //! The variables we want to save data
    std::set<std::string> _plotvariables;


    //! The filename suffix
    std::string _filename_suffix;

    
    //! Create the device
    void create_device(void);

    
    //! Create all the materials
    void create_materials(void);


    //! Create and setup the models
    void setup_models(void) throw (ModelErrorException);


};


//
// inline members
//


inline
Control::Control(const std::string& inputfile)
  : _inputfile(inputfile),
    _device(0),
    _database(0),
    _outputdir("."),
    _filename_suffix("")
{
}


inline
Device&
Control::get_device(void)
{
  return *_device;
}


inline
const std::set<std::string>&
Control::get_plotvariables(void) const
{
  return _plotvariables;
}


inline
const std::string&
Control::get_output_dir(void) const
{
  return _outputdir;
}


inline
const std::string&
Control::get_filename_suffix(void) const
{
  return _filename_suffix;
}


inline
void
Control::clear_filename_suffix(void)
{
  _filename_suffix = "";
}



inline
void
Control::set_filename_suffix(const std::string& suffix)
{
  if (suffix[0] != '_')
    _filename_suffix = "_" + suffix;
  else
    _filename_suffix = suffix;
}



inline
void
Control::append_to_filename_suffix(const std::string& suffix)
{
  if (suffix[0] != '_')
    _filename_suffix += "_" + suffix;
  else
    _filename_suffix += suffix;
}



#endif // _CONTROL_H_
