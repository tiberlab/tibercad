#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#include "SimulationInterface.h"
#include "SimulationOptions.h"
#include "InitFailedException.h"

#include <set>
#include <vector>
#include <string>

// HDF5 imports
#include <hdf5.h>
#include <hdf5_hl.h>


class TBDLLOCAL Interpolation : public SimulationInterface
{

  public:

    //! Destructor
    virtual ~Interpolation(void);

    //! Create a Interpolation object
    static Interpolation* create(const ModelOptions& options);

  protected:

    //! Constructor
    Interpolation(const ModelOptions& options);

    /*! \copydoc SimulationInterface::do_init() */
    virtual void do_init(void);

    virtual std::pair<double, double> get_value_and_derivative_secure(
              ID value_id, std::map<ID, double> params, ID dvar_id = INVALID_ID);

    /*! \copydoc SimulationInterface::do_solve() */
    virtual void do_solve(void);

    //! Abused here for the inheritance of variables
    virtual void do_print_info(void);

    //! Get global quantities from all sub-simulations
    virtual void get_solution_secure(const libMesh::Elem* elem, 
                                           std::map<ID, std::vector<double> >& values,
                                     const std::vector<libMesh::Point>& points);

    virtual void parse_options(void);

    virtual void do_setup_solution_variables(void);

  private:

   class model
   {
     public:
       void add_file(std::string filename);
       void add_model(std::string modelname);
       void add_principal_variable(std::string princ_var);
       void add_variables(const std::vector<std::string>& variables);

       //call read_data only after add methods have been called in sequence
       void read_data(void);

       const std::vector<std::string>& variables(void)
         { return _variables; };

       const std::vector<double>& princ_var_values(void)
         { return _princ_var_values; };

       const std::vector<double>& var_values(unsigned int princ_var_index, unsigned int variable_id);

       double data_value(unsigned int princ_var_index, std::vector<unsigned int>& var_indices);
       void set_log_data(bool flag=false);
       void set_log(const std::vector<std::string>& logvars);
       const bool log_data(void);
       const bool log_princ_var(void);
       const std::vector<bool>& log_var(void);

     private:
       std::string _filename;
       std::string _modelname;
       std::string _princ_var;
       std::vector<std::string> _variables;

       unsigned int _princ_var_size;
       std::vector<std::vector<unsigned int>> _variables_size;
       std::vector<unsigned int> _data_size;
       hid_t _princ_var_type;
       std::vector<std::vector<hid_t>> _variables_type;
       std::vector<hid_t> _data_type;

       std::vector<double> _princ_var_values;
       std::vector<std::vector<std::vector<double>>> _variables_values;
       std::vector<std::vector<double>> _data_values;

       bool _log_data;
       bool _log_princ_var;
       std::vector<bool> _log_var;


       //convert a n-dim array index to a 1-dim index in conformity with hdf5 internal indexing
       unsigned int _get_h5_index(const std::vector<unsigned int>& indices, const std::vector<unsigned int>& sizes);
   };

   std::vector<unsigned int> _find_range(double& value, const std::vector<double>& vec);

   std::map<ID,model*> _models;
   ID _princ_ID;

};

//
// inline methods
//


inline
Interpolation*
Interpolation::create(const ModelOptions& options)
{
  return new Interpolation(options);
}

inline
const std::vector<double>&
Interpolation::model::var_values(unsigned int princ_var_index, unsigned int variable_id)
{
  if (princ_var_index >= _princ_var_size )
    princ_var_index = 0;

  if (variable_id >= _variables.size() )
    variable_id = 0;

  return _variables_values[princ_var_index][variable_id];
}

inline
double
Interpolation::model::data_value(unsigned int princ_var_index, std::vector<unsigned int>& var_indices)
{
  if (princ_var_index >= _princ_var_size )
    princ_var_index = 0;

  var_indices.resize(_variables.size(), 0);

  for (unsigned int i = 0; i < var_indices.size(); i++)
  {
    if (var_indices[i] >= _variables_size[princ_var_index][i] )
      var_indices[i] = 0;
  }

  unsigned int index = _get_h5_index(var_indices, _variables_size[princ_var_index]);
  return _data_values[princ_var_index][index];
}

inline
void
Interpolation::model::set_log_data(bool flag)
{
  _log_data = flag;
}

inline
void
Interpolation::model::set_log(const std::vector<std::string>& logvars)
{
  std::set<std::string> vars;
  for (auto var : logvars)
    vars.insert(var);

  _log_princ_var = (vars.count(_princ_var)) ? true : false;

  _log_var.resize(_variables.size());
  for (unsigned int i = 0; i < _variables.size(); i++)
    _log_var[i] = (vars.count(_variables[i])) ? true : false;

}

inline
const bool
Interpolation::model::log_data(void)
{
  return _log_data;
}

inline
const bool
Interpolation::model::log_princ_var(void)
{
  return _log_princ_var;
}

inline
const std::vector<bool>&
Interpolation::model::log_var(void)
{
  return _log_var;
}

#endif // _INTERPOLATION_H_









