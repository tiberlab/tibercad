// $Id$

#include "tibercad/math/SignalGenerator.h"
#include "tibercad/base/Variable.h"
#include "tibercad/base/InitFailedException.h"

using namespace std;

SignalGenerator::SignalGenerator(const ModelOptions& options) :
  TiberModelObject(options),
  _input(0.0)
{
}

SignalGenerator::~SignalGenerator(void)
{
}


SignalGenerator*
SignalGenerator::create(const ModelOptions& options)
{
  string name = options.get_name();

  SignalGenerator* sg = create_from_library<SignalGenerator>(name, options);
  if (sg == nullptr)
  {
    sg = create_from_library<SignalGenerator>("signals/" + name, options);
  }

  if (sg == nullptr)
    throw InitFailedException("Cannot create signal class \"" + name + "\"");

  
  return(sg);
}


void
SignalGenerator::update_dependent_variables(void)
{
  do_update_dependent_variables();
}


void
SignalGenerator::init(void)
{
  // call this first, or the call to the initializer will break
  do_init();

  if (has_parameter("input"))
    get_parameter("input", _input, false,
        initializer(&SignalGenerator::update_dependent_variables));
  else
    VariableValue::check_and_register("$time", _input, this,
        initializer(&SignalGenerator::update_dependent_variables));

}
  
