// $Id$


#include "SimulationOptions.h"
#include "ModelOptions.h"
#include "Variable.h"


namespace {

  //! A private class to define the temperature as variable
  class _Temperature : Variable
  {

    public:

      _Temperature(void) {};

      void init(const std::string& s, double defaultval);

    protected:

      virtual void set_variable_value(double value, ID id = 0);

      virtual double get_variable_value(ID id = 0);
      
  };

  _Temperature _temperature;

  
  void
  _Temperature::set_variable_value(double value, ID id)
  {
    SimulationOptions::temperature = value;
  }

  double
  _Temperature::get_variable_value(ID id)
  {
    return SimulationOptions::temperature;
  }

  void
  _Temperature::init(const std::string& s, double defaultval)
  {
    SimulationOptions::temperature = check_and_register(s, defaultval);
  }

}




double
SimulationOptions::temperature = 300.0;

double&
SimulationOptions::temp = temperature;

double&
SimulationOptions::T = temperature;

bool
SimulationOptions::incomplete_ionization = false;

int
SimulationOptions::_verbose = 1;



int
SimulationOptions::verbose(void)
{
  return _verbose;
}


void
SimulationOptions::initialize(const ModelOptions& opts)
{
  std::string temp = opts.get_option("temperature", "300.0");
  _temperature.init(temp, 300.0);

  incomplete_ionization = opts.get_option("incomplete_ionization", true);

  _verbose = opts.get_option("verbose", 1);
}
    
