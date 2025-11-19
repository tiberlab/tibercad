// $Id$

#include "TriangleWave.h"
#include "tibercad/base/Variable.h"
#include "tibercad/base/InitFailedException.h"

//#include <fstream>

#include "tibercad/module/TiberModule.h"

using namespace std;

TriangleWave::TriangleWave(const ModelOptions& options) :
  SignalGenerator(options),
  _output(0.0),
  _output_str(""),
  _period(1e-3),
  _mean(0.0),
  _amplitude(0.5),
  _max_harmonic(10)
{
}


TriangleWave::~TriangleWave(void)
{
}


TriangleWave*
TriangleWave::create(const ModelOptions& options)
{
  return(new TriangleWave(options));
}


void
TriangleWave::do_init(void)
{
  _output_str = get_options().get_raw_option_string("output_variable");
  if (!VariableValue::check_string(_output_str))
    throw InitFailedException("The name \"" + _output_str +
        "\" given as output variable in Signal block is not a valid variable name");

  _output_str = VariableValue::check_and_register(_output_str, _output, this);

  get_parameter("period", _period);

  get_parameter("mean", _mean);

  get_parameter("amplitude", _amplitude);

  _max_harmonic = get_option("harmonics", _max_harmonic);

}

void
TriangleWave::do_update_dependent_variables(void)
{
  double tmp = 0.0;
  double sign = -1.0;

  for (unsigned int i = 0; i < _max_harmonic; ++i)
  {
    sign = -sign;
    double n = 2*i + 1;
    tmp += sign / (n * n) * sin(2*M_PI*n*get_input() / _period);
  }
  _output = _amplitude * tmp + _mean;

  VariableValue::set_variable_value(_output_str, _output);

  //ofstream of("test.dat", ios_base::app);
  //of << get_input() << " " << _output << endl;
}

