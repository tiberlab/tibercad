// $Id$

#ifndef _SIMULATIONOPTIONS_H_
#define _SIMULATIONOPTIONS_H_

//! Common options for set of simulations
class SimulationOptions
{

  public:

    //! The ambient temperature
    static double temperature;
    static double& temp;
    static double& T;

  private:
    
    SimulationOptions(void) {};

};

#endif //_SIMULATIONOPTIONS_H_
