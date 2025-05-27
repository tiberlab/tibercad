// $Id$

#ifndef _SIMULATIONOPTIONS_H_
#define _SIMULATIONOPTIONS_H_

#include "tiber_dll.h"
#include <string>

class ModelOptions;

//! Common options for set of simulations
class SimulationOptions
{

  public:


    //! The ambient temperature in K
    static double temperature;
    static double& temp;
    static double& T;

    //! Whether or not to consider incomplete ionization
    static bool incomplete_ionization;


    //! Tell verbosity
    static int verbose(void);


    //! Initialize the simulation options
    static void initialize(const ModelOptions& opts) TBDLLOCAL;

    //! A path for temporary data
    static std::string scratch_path;

  private:

    SimulationOptions(void) TBDLLOCAL {};
    
    //! Level of verbosity
    static int _verbose;

};



#endif //_SIMULATIONOPTIONS_H_
