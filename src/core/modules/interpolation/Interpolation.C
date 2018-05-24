#include "Interpolation.h"
#include "SimulationEnvironment.h"
#include "Constants.h"
#include "SolveFailedException.h"

// C++ includes
#include <fstream>
#include <math.h>

#include "TiberModule.h"

using namespace std;

Interpolation::Interpolation(const ModelOptions& options)
  : SimulationInterface(options)
{

}

void
Interpolation::do_init(void)
{
}

pair<double, double>
Interpolation::get_value_and_derivative_secure(
              ID value_id, map<ID, double> params, ID dvar_id)
{
  model* mod = _models[value_id];

  if (params.size() <= mod->variables().size() )
    return make_pair(0.0, 0.0);

  const vector<double>& princ_var = mod->princ_var_values();
  double princ_x = params[_princ_ID];

  params.erase(_princ_ID);
  if (dvar_id == _princ_ID)
    dvar_id = INVALID_ID;

  vector<unsigned int> princ_ids = _find_range(princ_x, princ_var);
  if (princ_ids[0] == princ_ids[1])
    princ_ids.resize(1);

  if (mod->log_princ_var())
    princ_x = log(princ_x);

  unsigned int dim = params.size();

  vector<double> f(princ_ids.size());
  vector<double> df(princ_ids.size());

  for (unsigned int p = 0; p < princ_ids.size(); p++)
  {
    vector<vector<double>> vars(dim, vector<double>(0));
    vector<vector<unsigned int>> range_id(dim, vector<unsigned int>(2));
    vector<vector<double>> range(dim, vector<double>(2));
    double vol = 1.0;

    vector<double> x(dim);
    vector<double> x0(dim), x1(dim);

    for (unsigned int i = 0; i < dim; i++)
    {
      x[i] = params[i];
      vars[i] = mod->var_values(princ_ids[p], i);
      range_id[i] = _find_range(x[i], vars[i]);

      if (mod->log_var()[i])
      {
        x[i] = log(x[i]);
        range[i][0] = log( vars[i][range_id[i][0]] );
        range[i][1] = log( vars[i][range_id[i][1]] );
      }
      else
      {
        range[i][0] = vars[i][range_id[i][0]];
        range[i][1] = vars[i][range_id[i][1]];
      }
      vol *= range[i][1]-range[i][0];

      x0[i] = x[i];
      x1[i] = x[i];

      if ( i == dvar_id )
      {
        x0[dvar_id] = range[dvar_id][0];
        x1[dvar_id] = range[dvar_id][1];
      }
    }


    unsigned int n_point = pow(2, dim);
    vector<double> vols(n_point, 1.0);
    vector<double> vols0(n_point, 1.0);
    vector<double> vols1(n_point, 1.0);
    vector<double> F(n_point);

    for (unsigned int i = 0; i < n_point; i++)
    {
      vector<unsigned int> F_ids(dim);
      for (unsigned int j = 0; j < dim; j++)
      {
        unsigned int k = (i>>j) & 1;
        vols[i] *= range[j][k] - x[j];
        F_ids[j] = range_id[j][k];

        if (dvar_id != INVALID_ID)
        {
          vols0[i] *= range[j][k] - x0[j];
          vols1[i] *= range[j][k] - x1[j];
        }
      }
      vols[i] = fabs(vols[i]);
      //cout<<"vol = "<<vol<<" vols["<<i<<"] = "<<vols[i]<<endl;

      if (dvar_id != INVALID_ID)
      {
        vols0[i] = fabs(vols0[i]);
        vols1[i] = fabs(vols1[i]);

        //cout<<"vol = "<<vol<<" vols0["<<i<<"] = "<<vols0[i]<<" vols1["<<i<<"] = "<<vols1[i]<<endl;
      }

      F[i] = (mod->log_data()) ? 
               log( mod->data_value(princ_ids[p], F_ids) ) : 
                 mod->data_value(princ_ids[p], F_ids);
    }

    double ft, f0, f1, dft = 0.0;
    for (unsigned int i = 0; i < n_point; i++)
    {
      unsigned int j = n_point - 1 - i;
      ft += vols[i]*F[j];

      if (dvar_id != INVALID_ID)
      {
        f0 += vols0[i]*F[j];
        f1 += vols1[i]*F[j];
      }
    }

    ft /= vol;

    if (dvar_id != INVALID_ID)
    {
      f0 /= vol;
      f1 /= vol;
      dft = (f1-f0)/(x1[dvar_id]-x0[dvar_id]);

      if (mod->log_var()[dvar_id])
        dft /= exp(x[dvar_id]);
    }

    f[p] = ft;
    df[p] = dft;

  } //principal variable loop

  if (f.size() > 1)
  {
    double p0 = (mod->log_princ_var()) ? log(princ_var[princ_ids[0]]) : princ_var[princ_ids[0]];
    double p1 = (mod->log_princ_var()) ? log(princ_var[princ_ids[1]]) : princ_var[princ_ids[1]];

    f[0] = (f[0]*(p1 - princ_x) + f[1]*(princ_x - p0))/(p1-p0);
    df[0] = (df[0]*(p1 - princ_x) + df[1]*(princ_x - p0))/(p1-p0);
  }

  if (mod->log_data())
  {
    f[0] = exp(f[0]);
    df[0] = f[0]*df[0];
  }

  pair<double, double> value = make_pair(f[0], df[0]);
  return value;
}

void
Interpolation::do_solve(void)
{
}

void
Interpolation::do_print_info(void)
{
}

void
Interpolation::get_solution_secure(const libMesh::Elem* elem, 
                                           map<ID, vector<double> >& values,
                                     const vector<libMesh::Point>& points)
{
}

void
Interpolation::parse_options(void)
{
}

void
Interpolation::do_setup_solution_variables(void)
{
  auto itm (get_options().submodels_begin("model"));
  auto end_itm (get_options().submodels_end("model"));

  for ( ; itm != end_itm; ++itm)
  {
    string name = (itm->second).get_option("name", "");
    string file = (itm->second).get_option("file_name", "");
    string princvar = (itm->second).get_option("principal_variable", "");
    vector<string> variables;
    (itm->second).get_option("variables", variables);

    if (name == "")
      throw InitFailedException("Interpolation: model name must be specified");
    if (file == "")
      throw InitFailedException("Interpolation: file name for model '" + name + "'must be specified");
    if (princvar == "")
      throw InitFailedException("Interpolation: principal variable for model '" + name + "'must be specified");
    if (variables.size() == 0)
      throw InitFailedException("Interpolation: variables for model '" + name + "'must be specified");

    model* mod = new model;

    mod->add_file(file);
    mod->add_model(name);
    mod->add_principal_variable(princvar);
    mod->add_variables(variables);

    mod->read_data();

    bool logdata = (itm->second).get_option("log_data", false);
    vector<string> logvars;
    (itm->second).get_option("log_variables", logvars);
    mod->set_log_data(logdata);
    mod->set_log(logvars);

    add_value(name);
    ID id = get_value_id(name);
    _models.insert( make_pair(id,mod) );

    //parameter ids will be the same as inside mod
    for (auto var : mod->variables())
      add_parameter(var, id);

    //last id will be principal variable
    add_parameter(princvar, id);
    _princ_ID = get_param_id(princvar, id);

  }

}

Interpolation::~Interpolation(void)
{
}



vector<unsigned int>
Interpolation::_find_range(double& value, const vector<double>& vec)
{
  vector<unsigned int> ids(2);

  if (vec.size() == 1)
  {
    value = vec[0];
    ids[0] = 0;
    ids[1] = 0;
    return ids;
  }

  unsigned int idx1, idx2;

  if (value < vec.front())  // check if value is in the whole range, otherwise don't even try to search
  {
    value = vec.front();
    idx1 = 0;
    idx2 = 1;
  }
  else if (value > vec.back())
  {
    value = vec.back();
    idx1 = vec.size() - 2;
    idx2 = vec.size() - 1;
  }
  else  // find the range
  {
    vector<double> _tmp = vec;
    idx1 = 0;
    idx2 = _tmp.size() - 1;
    while (_tmp.size() > 2) {
      size_t half_size = _tmp.size() / 2;
      vector<double> _lo(_tmp.begin(), _tmp.begin() + half_size);
      vector<double> _hi(_tmp.begin() + half_size, _tmp.end());

      if ( value >= _hi[0])
      {
        idx1 += half_size;
        _tmp = _hi;
      }
      else
      {
        idx2 = idx1 + half_size;
        _lo.push_back(_hi[0]);
        _tmp = _lo;
      }
    } // end while
  } // end else

  ids[0] = idx1;
  ids[1] = idx2;

  return ids;
}




void
Interpolation::model::add_file(std::string filename)
{
  hid_t file_id;
  file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file_id < 0)
    throw InitFailedException("'" + filename + "' is not a valid hdf5 file");

  _filename = filename;

  H5Fclose(file_id);
}


void
Interpolation::model::add_model(std::string modelname)
{
  hid_t file_id, model_id;
  file_id = H5Fopen(_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

  model_id = H5Gopen2(file_id, modelname.c_str(), H5P_DEFAULT);

  if (model_id < 0)
    throw InitFailedException("Model '" + modelname + "' not found in file '" + _filename + "'");

  _modelname = modelname;

  H5Gclose(model_id);
  H5Fclose(file_id);
}


void
Interpolation::model::add_principal_variable(std::string princ_var)
{
  hid_t file_id, model_id, princvar_id, data_type, native_type;
  hsize_t* dims;
  int rank;

  file_id = H5Fopen(_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  model_id = H5Gopen2(file_id, _modelname.c_str(), H5P_DEFAULT);

  princvar_id = H5Dopen2(model_id, princ_var.c_str(), H5P_DEFAULT);
  // Get the data type of the dataset
  data_type = H5Dget_type(princvar_id);
  native_type = H5Tget_native_type(data_type, H5T_DIR_ASCEND);

  if (princvar_id < 0)
    throw InitFailedException("Principal variable '" + princ_var + "' not found\nfor model '" + _modelname + "' in file '" + _filename + "'");

  _princ_var = princ_var;
  _princ_var_type = native_type;

  //get dimensionality
  H5LTget_dataset_ndims(model_id, _princ_var.c_str(), &rank);

  if (rank > 1)
    throw InitFailedException("Principal variable '" + _princ_var + "' rank is not 1.\nError in model '" + _modelname + "' in file '" + _filename + "'");

  //get size
  dims = (hsize_t*)malloc(sizeof(hsize_t) * rank);
  H5LTget_dataset_info(model_id, _princ_var.c_str(), dims, NULL, NULL);
  _princ_var_size = dims[0];
  free(dims);

  H5Dclose(princvar_id);

  for (unsigned int i = 0; i < _princ_var_size; i++)
  {
    std::string gname = _princ_var + "_" + std::to_string(i);
    princvar_id = H5Gopen2(model_id, gname.c_str(), H5P_DEFAULT);

    if (princvar_id < 0)
      throw InitFailedException("Data groups of model '" + _modelname + "' in file '" + _filename + "'\n are not consistent with principal variable '" + _princ_var + "' vector size");

    H5Gclose(princvar_id);
  }

  H5Gclose(model_id);
  H5Fclose(file_id);
}

void
Interpolation::model::add_variables(const std::vector<std::string>& variables)
{
  hid_t file_id, model_id, princvar_id, var_id, attr_id, data_id, data_type, native_type;
  hsize_t* dims;
  int rank;
  size_t  type_size;
  H5T_class_t type_class;

  file_id = H5Fopen(_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  model_id = H5Gopen2(file_id, _modelname.c_str(), H5P_DEFAULT);

  _variables_size.resize(0);
  _variables_type.resize(0);
  _data_size.resize(0);
  _data_type.resize(0);
  _variables_size.resize( _princ_var_size, vector<unsigned int>(variables.size()) );
  _variables_type.resize( _princ_var_size, vector<hid_t>(variables.size()) );
  _data_size.resize(_princ_var_size);
  _data_type.resize(_princ_var_size);


  // loop over principal variable groups
  for (unsigned int i = 0; i < _princ_var_size; i++)
  {
     _variables.resize(0);
     _variables.resize(variables.size(), "");

     // open the group associated to each principal variable value
     string gname = _princ_var + "_" + to_string(i);
     princvar_id = H5Gopen2(model_id, gname.c_str(), H5P_DEFAULT);

     unsigned int data_size = 1;

     // loop over variables of each group
     for (unsigned int j = 0; j < variables.size(); j++)
     {
       hsize_t* vdims;

       // open the variable vector to check it exists
       var_id = H5Dopen2(princvar_id, variables[j].c_str(), H5P_DEFAULT);

       if (var_id < 0)
         throw InitFailedException("Variable '" + variables[j] + "' not found\nfor model '" + _modelname + "' in file '" + _filename + "'");

       // Get the data type of the dataset
       data_type = H5Dget_type(var_id);
       native_type = H5Tget_native_type(data_type, H5T_DIR_ASCEND);


       //get dimensionality of variable vector, rank should be 1
       H5LTget_dataset_ndims(princvar_id, variables[j].c_str(), &rank);

       if (rank > 1)
         throw InitFailedException("Variable '" + variables[j] + "' rank is not 1.\nError in model '" + _modelname + "' in file '" + _filename + "'");

       //get size of each variable vector and multiplicate to get the total size of model 'data' values
       vdims = (hsize_t*)malloc(sizeof(hsize_t)*rank);
       H5LTget_dataset_info(princvar_id, variables[j].c_str(), vdims, NULL, NULL);

       data_size *= vdims[0];  //at the end data_size will be compared with the size of 'data' dataset in the file
       unsigned int var_size = vdims[0];
       hid_t var_type = native_type;

       //get variable ids checking they are unique and not greater than the number of variables
       attr_id = H5Aopen(var_id, "id", H5P_DEFAULT);

       if (attr_id < 0)
         throw InitFailedException("Attribute 'id' for variable '" + variables[j] + "' not found.\nError in model '" + _modelname + "' in file '" + _filename + "'");

       //get attribute dimensionality and check it is 1
       H5LTget_attribute_ndims(princvar_id, variables[j].c_str(), "id", &rank );

       if (rank > 1)
         throw InitFailedException("Attribute 'id' for variable '" + variables[j] + "' should have rank 1.\nError in model '" + _modelname + "' in file '" + _filename + "'");

       // check attribute size is 1
       vdims = (hsize_t*)realloc(vdims, sizeof(hsize_t)*rank);
       H5LTget_attribute_info(princvar_id, variables[j].c_str(), "id", vdims, &type_class, &type_size);

       if (vdims[0] > 1)
         throw InitFailedException("Attribute 'id' for variable '" + variables[j] + "' should have size 1.\nError in model '" + _modelname + "' in file '" + _filename + "'");

       // read the attribute "id" and check it is unique and not greater than the number of variables
       // Get the data type of the attribute
       data_type = H5Aget_type(attr_id);
       native_type = H5Tget_native_type(data_type, H5T_DIR_ASCEND);
       // Allocate and read data
       int* id = (int *)malloc(vdims[0] * sizeof(int));
       H5Aread(attr_id, native_type, id );
       //H5Aread(attr_id, H5T_NATIVE_INT, id );

       if (id[0] > (variables.size() - 1))
         throw InitFailedException("Attribute 'id' for variable '" + variables[j] + "' exceeds the number of variables.\nError in model '" + _modelname + "' in file '" + _filename + "'");

       if (_variables[id[0]] != "")
         throw InitFailedException("Attribute 'id' for variable '" + variables[j] + "' is not unique.\nError in model '" + _modelname + "' in file '" + _filename + "'");

       _variables[id[0]] = variables[j];
       _variables_size[i][id[0]] = var_size;
       _variables_type[i][id[0]] = var_type;

       H5Aclose(attr_id);
       H5Dclose(var_id);
       free(vdims);
       free(id);
     }

     //open 'data'
     data_id = H5Dopen2(princvar_id, "data", H5P_DEFAULT);

     if (data_id < 0)
         throw InitFailedException("'data' not found for model '" + _modelname + "' in file '" + _filename + "'");

     // Get the data type of the dataset
     data_type = H5Dget_type(data_id);
     native_type = H5Tget_native_type(data_type, H5T_DIR_ASCEND);

     //get dimensionality of variable vector
     H5LTget_dataset_ndims(princvar_id, "data", &rank);

     if (rank != variables.size())
       throw InitFailedException("'data' rank not consistent with the number of variables.\nError in model '" + _modelname + "' in file '" + _filename + "'");

     dims = (hsize_t*)malloc(sizeof(hsize_t)*rank);
     H5LTget_dataset_info(princvar_id, "data", dims, NULL, NULL);
     unsigned int tot_dim = 1;
     for (unsigned int j = 0; j < rank; j++)
       tot_dim *= dims[j];

     free(dims);

     if (data_size != tot_dim)
       throw InitFailedException("'data' size not consistent with variables sizes.\nError in model '" + _modelname + "' in file '" + _filename + "'");

     _data_type[i] = native_type;
     _data_size[i] = tot_dim;

     H5Dclose(data_id);
     H5Gclose(princvar_id);
  }
  H5Gclose(model_id);
  H5Fclose(file_id);


}


void
Interpolation::model::read_data(void)
{
  _princ_var_values.resize(0);
  _princ_var_values.resize(_princ_var_size);
  _variables_values.resize(0);
  _variables_values.resize( _princ_var_size, vector<vector<double>>(0) );
  _data_values.resize(0);
  _data_values.resize(_princ_var_size, vector<double>(0) );

  hid_t file_id, model_id, princvar_id, var_id, attr_id, data_id;

  file_id = H5Fopen(_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  model_id = H5Gopen2(file_id, _modelname.c_str(), H5P_DEFAULT);

  // read the values of the principal variable
  princvar_id = H5Dopen2(model_id, _princ_var.c_str(), H5P_DEFAULT);
  H5Dread(princvar_id, _princ_var_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, _princ_var_values.data());
  H5Dclose(princvar_id);

  // loop over principal variable groups
  for (unsigned int i = 0; i < _princ_var_size; i++)
  {
    // open the group associated to each principal variable value
    string gname = _princ_var + "_" + to_string(i);
    princvar_id = H5Gopen2(model_id, gname.c_str(), H5P_DEFAULT);

    _data_values[i].resize(_data_size[i]);

    //read 'data'
    data_id = H5Dopen2(princvar_id, "data", H5P_DEFAULT);
    H5Dread(data_id, _data_type[i], H5S_ALL, H5S_ALL, H5P_DEFAULT, _data_values[i].data());
    H5Dclose(data_id);

    _variables_values[i].resize(_variables.size(), vector<double>(0) );

    // loop over variables of each group
    for (unsigned int j = 0; j < _variables.size(); j++)
    {
      _variables_values[i][j].resize(_variables_size[i][j]);

      //read variable vector
      var_id = H5Dopen2(princvar_id, _variables[j].c_str(), H5P_DEFAULT);
      H5Dread(var_id, _variables_type[i][j], H5S_ALL, H5S_ALL, H5P_DEFAULT, _variables_values[i][j].data());
      H5Dclose(var_id);
    }

    H5Gclose(princvar_id);
  }

  H5Gclose(model_id);
  H5Fclose(file_id);



/* DATA OUTPUT
  cout<<"PRINCIPAL VARIABLE: "<<_princ_var<<endl;
  cout<<"values: "<<endl;
  for (auto val : _princ_var_values)
    cout<<val<<endl;
  cout<<endl;

  for (unsigned int i = 0; i < _princ_var_size; i++)
  {
    cout<<"VARIABLES"<<endl;
    cout<<_princ_var<<"_"<<i<<endl;
    for (unsigned int j = 0; j < _variables[i].size(); j++)
    {
      cout<<"variable["<<j<<"]: "<<_variables[i][j]<<endl;
      for (unsigned int k = 0; k < _variables_size[i][j]; k++)
        cout<<"["<<k<<"] = "<<_variables_values[i][j][k]<<endl;
    }
    cout<<"DATA"<<endl;
    for (unsigned int j = 0; j < _data_size[i]; j++)
    cout<<"["<<j<<"] = "<<_data_values[i][j]<<endl;
  }
*/
}


unsigned int 
Interpolation::model::_get_h5_index(const vector<unsigned int>& indices, const vector<unsigned int>& sizes)
{
  unsigned int n = indices.size();
  unsigned int index = indices[n-1];

  for (size_t i = 0; i <= n-2; i++)
  {
    unsigned int N = 1;
    for (size_t j = i+1; j <= n-1; j++)
      N *= sizes[j];

    index += indices[i]*N;
  }

  return index;
}







