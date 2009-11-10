// $Id$

#ifndef _SIMULATIONOPTIONS_H_
#define _SIMULATIONOPTIONS_H_


class ModelOptions;

//! Common options for set of simulations
class SimulationOptions
{

  public:

    SimulationOptions(void) {};


    //! The ambient temperature
    static double temperature;
    static double& temp;
    static double& T;


    //! Wheter or not to consider incomplete ionization
    static bool incomplete_ionization;


    //! Tell verbosity
    static int verbose(void);


    //! Initialize the simulation options
    static void initialize(const ModelOptions& opts);


  private:
    
    //! Level of verbosity
    static int _verbose;

};



#endif //_SIMULATIONOPTIONS_H_
