// $Id$

#include "SignalGenerator.h"
#include "Variable.h"

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

  
}


void
SignalGenerator::update_dependent_variables(void)
{
  do_update_dependent_variables();
}


void
SignalGenerator::init(void)
{
  if (has_parameter("input"))
    get_parameter("input", _input);
  else
    VariableValue::check_and_register("$time", _input, this,
        initializer(&SignalGenerator::update_dependent_variables));

  do_init();
}
  
