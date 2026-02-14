/*  
 * This file is part of the tiberCAD module tmm.
 *
 * tiberCAD modules are licensed under the GNU General Public License v3.
 *
 * tiberCAD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * tiberCAD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tiberCAD. If not, see <http://www.gnu.org/licenses/>.
 */

/*!
 * \file Tmm.C
 * \brief tiberCAD tmm module implementation.
 *
 * \note This file is part of module tmm.
 */


#include "Tmm.h"
#include "tibercad/solver/TiberLinearSystem.h"
#include "TmmBulkModel.h"
#include "TmmBoundaryModel.h"
#include "tibercad/module/SimulationEnvironment.h"
#include "tibercad/base/Device.h"
#include "tibercad/io/Messages.h"
#include "tibercad/io/Database.h"
#include <filesystem>
#include <cctype>
#include <cmath>
#include "matrix2by2.h"
#include <algorithm>
#include "libmesh/dof_map.h"




// This is needed in order to create the shared module library
#include "tibercad/module/TiberModule.h"


using namespace libMesh;



// ,_incident_angle({0.0})
Tmm::Tmm(const ModelOptions& options) :
  SimulationInterface(options)
{
}

std::complex<double> Tmm::cmlx_sqrt(std::complex<double> in)
{
  double real_part = real(in);
  double imag_part = imag(in);

  if (imag_part == 0)
    imag_part = 0.0;
  if (real_part == 0)
    real_part = 0.0;
  std::complex<double> out(real_part, imag_part);
  out = sqrt(out);


  return(out);

}


bool Tmm::cmp_string(const std::string& a, const std::string& b)
{
  unsigned int sz = a.size();
  if (b.size() != sz)
    return false;
  for (unsigned int i = 0; i < sz; ++i)
    if (tolower(a[i]) != tolower(b[i]))
      return false;
  return true;
}



matrix2by2 Tmm::get_M(double n_real, double n_imag, double lenght, double lambda, double kr, double phase)
{
  matrix2by2 new_Matrix_2by2;
  complex<double> j(0, 1);
  complex<double> ki((2 * M_PI * n_real) / lambda, (2 * M_PI * n_imag) / lambda);
  complex<double> kz;
  //complex<double> memory;
  kz = ki * cmlx_sqrt(complex<double>(1.0 - (pow(kr, 2) / pow(ki, 2))));
  //std::cout<<"kr is :"<<kr<<"   ki :"<<ki<< "   kz :"<<kz<<"    "<<exp(-j * kz)<<std::endl;
  kz *= lenght;
  complex<double> ps(0, phase);

  new_Matrix_2by2.set(0, exp((-j * kz) + ps));
  new_Matrix_2by2.set(1, 0);
  new_Matrix_2by2.set(2, 0);
  new_Matrix_2by2.set(3, exp((j * kz) - ps));


  return(new_Matrix_2by2);
}

matrix2by2 Tmm::get_D(double n1_real, double n1_imag, double n2_real, double n2_imag, double kr, double lambda, double mode)
{
  matrix2by2 new_Matrix_2by2;
  complex<double> j(0, 1);
  complex<double> Ni1(n1_real, n1_imag);
  complex<double> Ni2(n2_real, n2_imag);
  complex<double> ki1((2 * M_PI * n1_real) / lambda, (2 * M_PI * n1_imag) / lambda);
  complex<double> ki2((2 * M_PI * n2_real) / lambda, (2 * M_PI * n2_imag) / lambda);
  ki1 = ki1 * cmlx_sqrt(complex<double>(1.0 - (pow(kr, 2) / pow(ki1, 2))));
  ki2 = ki2 * cmlx_sqrt(complex<double>(1.0 - (pow(kr, 2) / pow(ki2, 2))));
  complex<double> r, t;


  if (mode == 0)  //TE Mode
  {
    r = (ki1 - ki2) / (ki1 + ki2);
    t = (2.0 * ki1) / (ki1 + ki2);
  }
  else           //TM Mode
  {
    r = ((-ki1 * Ni2 * Ni2) + (ki2 * Ni1 * Ni1)) / ((ki1 * Ni2 * Ni2) + (ki2 * Ni1 * Ni1));
    t = (2 * ki1 * Ni1 * Ni2) / ((ki1 * Ni2 * Ni2) + (ki2 * Ni1 * Ni1));
  }

  new_Matrix_2by2.set(0, 1.0 / (t));
  new_Matrix_2by2.set(1, (r) / (t));
  new_Matrix_2by2.set(2, (r) / (t));
  new_Matrix_2by2.set(3, 1.0 / (t));

  return(new_Matrix_2by2);
}

matrix2by2 Tmm::Determinal_Matrix(matrix2by2 MAT)
{
  complex<double> ratio;
  matrix2by2 MAT_DET;
  ratio = 1 / (MAT.get(0) * MAT.get(3) - MAT.get(1) * MAT.get(2));
  MAT_DET.set(0, ratio * MAT.get(3));
  MAT_DET.set(1, -ratio * MAT.get(1));
  MAT_DET.set(2, -ratio * MAT.get(2));
  MAT_DET.set(3, ratio * MAT.get(0));
  return (MAT_DET);

}

void Tmm::reset_global_variables()
{
  _Transmission.clear();
  _Reflection.clear();
  _Absorption.clear();
  _AVT = 0;

  _Electric_Field_External.clear();
  _Electric_Field_Internal.clear();
  _Generation_regions.clear();
  _Generation_regions_internal.clear();
  _regions_name.clear();

  _Intensity.clear();
  _External_Source_ElectricField.clear();
  _Wavelength.clear();
  _Generation_rate.clear();
  _Poynting_external.clear();
  _Energy_loss_external.clear();

  _Internal_Wavelength.clear();
  _Internal_Source_ElectricField.clear();
  _Internal_Poynting.clear();
  _Internal_Power.clear();
  _Internal_Absorption.clear();
  _Abs.clear();
  _Internal_Intensity.clear();

  _angle.clear();
  _Poynting_front.clear();
  _Poynting_back.clear();

  _solutions.clear();

  _green_vector.clear();
}

void Tmm::dipole_source(double& A_P, double& A_N, double Mode, double Oraintation, double lambda, double cos_phi_inter)
{

  double sin_phi_inter;
  sin_phi_inter = real(cmlx_sqrt(complex<double>(1.0 - pow(cos_phi_inter, 2))));


  if (Mode == 0)    //TE Mode
  {
    if (Oraintation == 0)  // Horizontal 
    {
      A_P = -sqrt(3 / (16 * M_PI));
      A_N = sqrt(3 / (16 * M_PI));
    }
    else   //Vertical
    {
      A_P = -sqrt(3 / (16 * M_PI)) * sin_phi_inter;
      A_N = sqrt(3 / (16 * M_PI)) * sin_phi_inter;
    }
  }
  else     //TM Mode
  {
    if (Oraintation == 0) // Horizontal 
    {
      A_P = -sqrt(3 / (16 * M_PI)) * cos_phi_inter;
      A_N = sqrt(3 / (16 * M_PI)) * cos_phi_inter;
    }
    else   //Vertical
    {
      A_P = -sqrt(3 / (16 * M_PI)) * sin_phi_inter;
      A_N = sqrt(3 / (16 * M_PI)) * sin_phi_inter;
    }

  }
}

void Tmm::solving_internal_source(vector<complex<double>>& E_int_f, vector<complex<double>>& E_int_b, const vector<double>& n_real,
  const vector<double>& n_imag, const vector<double>& l, double lambda, double dipole_loc, double kr, double Mode, double A_P, double phase)
{
  matrix2by2 Es(A_P, 0, -A_P, 0);
  matrix2by2 DD(1, 0, 0, 1);
  matrix2by2 MM(0, 0, 0, 0);
  matrix2by2 TT_load(0, 0, 0, 0);
  matrix2by2 T_RIGHT(1, 0, 0, 1);
  matrix2by2 T_LEFT(1, 0, 0, 1);
  matrix2by2 T_TOTAL(1, 0, 0, 1);
  matrix2by2 AA(0, 0, 0, 1);


  size_t mesh_size = n_real.size();
  //*****************************loop for right side of the dipole******************************
  MM = get_M(n_real[mesh_size - 1], n_imag[mesh_size - 1], l[mesh_size - 1], lambda, kr, phase);

  T_RIGHT = MM;
  for (size_t k = mesh_size - 1; k > dipole_loc; --k)
  {
    DD = get_D(n_real[k - 1], n_imag[k - 1], n_real[k], n_imag[k], kr, lambda, Mode);
    MM = get_M(n_real[k - 1], n_imag[k - 1], l[k - 1], lambda, kr, phase);
    TT_load = DD * T_RIGHT;
    T_RIGHT = MM * TT_load;
  }
  MM = get_M(n_real[dipole_loc - 1], n_imag[dipole_loc - 1], l[dipole_loc - 1] / 2, lambda, kr, phase);
  TT_load = MM * T_RIGHT;
  T_RIGHT = TT_load;

  //*****************************loop for left side of the dipole******************************
  MM = get_M(n_real[dipole_loc - 1], n_imag[dipole_loc - 1], l[dipole_loc - 1] / 2, lambda, kr, phase);
  T_LEFT = MM;
  for (size_t k = dipole_loc - 2; k > 0; --k)
  {
    DD = get_D(n_real[k], n_imag[k], n_real[k + 1], n_imag[k + 1], kr, lambda, Mode);
    MM = get_M(n_real[k], n_imag[k], l[k], lambda, kr, phase);
    TT_load = DD * T_LEFT;
    T_LEFT = MM * TT_load;
  }
  MM = get_M(n_real[0], n_imag[0], l[0], lambda, kr, phase);
  TT_load = MM * T_LEFT;
  T_LEFT = TT_load;
  // T_LEFT.print();
  //*****************************Solving a pair of linear equations to get E.F. in two ends******************************
  T_LEFT.inv();
  matrix2by2 Sol(1, 0, 0, 1);
  Sol.set(0, T_RIGHT.get(0));
  Sol.set(1, -T_LEFT.get(1));
  Sol.set(2, T_RIGHT.get(2));
  Sol.set(3, -T_LEFT.get(3));
  Sol.inv();
  TT_load = Sol * Es;                 // Electric Field in two ends
  // std::cout<<"Out side E.F is:"<<std::endl;
  //  TT_load.print();
  //*****************************Defining vectors for Electric Field******************************
  vector<complex<double>> E_int_f_r(mesh_size);
  vector<complex<double>> E_int_b_r(mesh_size);
  vector<complex<double>> E_int_f_l(mesh_size);
  vector<complex<double>> E_int_b_l(mesh_size);
  //*****************************Assigning Electric field of two ends of the device******************************
  E_int_f_r[mesh_size - 1] = TT_load.get(0);
  E_int_b_r[mesh_size - 1] = 0;
  E_int_f_l[0] = 0;
  E_int_b_l[0] = TT_load.get(2);
  //*********using last sections results to calculate electric field for the entire device******************
  matrix2by2 E_N(1, 0, 0, 0);
  matrix2by2 E_I(0, 0, 0, 0);
  E_N.set(0, E_int_f_r[mesh_size - 1]);
  E_N.set(2, E_int_b_r[mesh_size - 1]);

  T_RIGHT.unit_matrix();
  DD.unit_matrix();
  //*****************************loop for right side TMM******************************
  for (size_t k = mesh_size; k > dipole_loc; --k)
  {
    if (k < mesh_size)
      DD = get_D(n_real[k - 1], n_imag[k - 1], n_real[k], n_imag[k], kr, lambda, Mode);
    MM = get_M(n_real[k - 1], n_imag[k - 1], l[k - 1], lambda, kr, phase);
    TT_load = DD * T_RIGHT;
    T_RIGHT = MM * TT_load;
    E_I = T_RIGHT * E_N;
    E_int_f_r[k - 2] = (E_I.get(0));
    E_int_b_r[k - 2] = (E_I.get(2));
  }

  E_N.set(0, E_int_b_l[0]);
  E_N.set(2, E_int_f_l[0]);


  T_LEFT.unit_matrix();
  DD.unit_matrix();
  //*****************************loop for left side TMM******************************
  for (size_t k = 0; k < dipole_loc - 1; ++k)
  {
    if (k > 0)
      DD = get_D(n_real[k - 1], n_imag[k - 1], n_real[k], n_imag[k], kr, lambda, Mode);
    MM = get_M(n_real[k], n_imag[k], l[k], lambda, kr, phase);
    DD.inv();
    TT_load = DD * T_LEFT;
    T_LEFT = MM * TT_load;
    E_I = T_LEFT * E_N;
    E_int_f_l[k + 1] = (E_I.get(2));
    E_int_b_l[k + 1] = (E_I.get(0));
  }

  //*****************************loops to to sum up electric fields******************************

  for (size_t nm = 0; nm < dipole_loc;nm++)
  {
    E_int_f[nm] = E_int_f_l[nm + 1];
    E_int_b[nm] = E_int_b_l[nm + 1];
  }

  for (size_t nm = E_int_f_r.size() - 1; nm > dipole_loc - 1;nm--)
  {
    E_int_f[nm] = E_int_f_r[nm - 1];
    E_int_b[nm] = E_int_b_r[nm - 1];
  }
  E_int_f[dipole_loc - 1] = (E_int_f_l[dipole_loc - 1] + E_int_f_r[dipole_loc - 1]) / 2;
  E_int_b[dipole_loc - 1] = (E_int_b_l[dipole_loc - 1] + E_int_b_r[dipole_loc - 1]) / 2;
}



vector<double> Tmm::linear_interpolation1(vector<double> xData, vector<double> yData, vector<double> x_interp)
{
  int size = xData.size();
  vector<double> y_interp;

  for (double x : x_interp)
  {
    int i = 0;
    if (x >= xData[size - 2])
    {
      i = size - 2;
    }
    else
    {
      while (x > xData[i + 1]) i++;
    }
    double xL = xData[i], yL = yData[i], xR = xData[i + 1], yR = yData[i + 1];
    double delta = (yR - yL) / (xR - xL);
    double y = yL + delta * (x - xL);
    y_interp.push_back(y);
  }
  return y_interp;
}



Tmm::~Tmm(void)
{
  string outdir = get_output_directory();
  std::filesystem::remove(outdir + "/green_matrix.dat");
}


Tmm*
Tmm::create(const ModelOptions& options)
{
  return new Tmm(options);
}



void
Tmm::do_init(void)
{
  parse_options();

  // _env = &get_environment();
  // _device = &(_env->get_device());

 // create a linear equation system 
  create_equation_system("linear");

  // get the reference to it
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  // add variables and attach the assemble function
  system.add_variable("G", libMeshEnums::CONSTANT, MONOMIAL, &get_region_ids());



  system.init();
}


void
Tmm::parse_options(void)
{
  // read wavelengths from input

  _incident_angle = get_option("incident_angle", 0);

  _reflectivity = get_option("back_reflectivity", 0);

  _dipole_loc = get_option("dipole_loc", 0);


  _coh_mod = get_option("coh_mode", 0);



  _polarization = get_option("polarization", "");
  _orientation = get_option("orientation", "");

  if (_polarization.empty())
    _polarization_vec.push_back(0);

  if (cmp_string(_polarization, "TEM"))
  {
    _polarization_vec.push_back(0);
    _polarization_vec.push_back(1);
  }
  else
  {
    if (cmp_string(_polarization, "TE"))
      _polarization_vec.push_back(0);
    else
      if (cmp_string(_polarization, "TM"))
        _polarization_vec.push_back(1);
      else
        Messages::warning("Unknown polarization.");
  }



  if (_orientation.empty())
    _oraintation_vec.push_back(0);

  if (cmp_string(_orientation, "VH"))
  {
    _oraintation_vec.push_back(0);
    _oraintation_vec.push_back(1);
  }
  else
  {
    if (cmp_string(_orientation, "H"))
      _oraintation_vec.push_back(0);
    else
      if (cmp_string(_orientation, "V"))
        _oraintation_vec.push_back(1);
      else
        Messages::warning("Unknown orientation.");
  }


  _steps = get_option("wave_number_steps", 0);
  get_option("wave_number_ratio", _ratio);





  get_option("dipole_power", _dipole_power);
  get_option("dipole_coordinate", _dipole_coordinate);

  if (_dipole_coordinate.size() != _dipole_power.size())
    if (_dipole_power.size() == 1)
    {
      _dipole_power.resize(_dipole_coordinate.size());
      for (int nm = 1; nm < _dipole_coordinate.size(); nm++)
        _dipole_power[nm] = _dipole_power[0];
    }
    else
      Messages::warning("The size of dipole_power and dipole_coordinate vectors don't match!");

  get_option("wavelengths", _wavelength_vector);
  if (_wavelength_vector.empty())
  {
    _up_lambda = get_option("wavelength_uper_lim", 0);
    if (_up_lambda == 0)
    {
      Messages::warning("You did not provide any upper wavelength limit for TMM.");
    }
    _down_lambda = get_option("wavelength_lower_lim", 0);
    if (_down_lambda == 0)
    {
      Messages::warning("You did not provide any lower wavelength limit for TMM.");
    }

    _wavelength_steps = get_option("wavelength_steps", 1);
  }


  Database db;
  ifstream is;
  db.set_material("Sun1p5am", get_option("illumination_spectrum", ""));
  is.open(db.get_data_file().c_str());
  if (!is.fail() || is.good())
  {
    //throw InitFailedException("Cannot read spectrum "
      //  "from file " + db.get_data_file());

    size_t i = 0;
    const size_t buf_len = 256;
    char buf[buf_len];

    while (is.good())
    {
      if (i == _lambda.size())
      {
        size_t n_new = _lambda.size() + 100;
        _lambda.reserve(n_new);
        _spectrum.reserve(n_new);
      }

      is.getline(buf, buf_len);
      if (buf[0] != '#')
      {
        istringstream in(buf);

        double l, s;
        if (in >> l >> s)
        {
          _lambda.push_back(l);
          // conversion from nm^-1 to J
          _spectrum.push_back(s);
          i++;
        }
      }
    }
    is.close();
  }
  else
    Messages::warning("Couldn't find External source spectrum");

  Database datab;
  ifstream ifs;
  datab.set_material("Eye_sens", get_option("Eye_sens", ""));
  ifs.open(datab.get_data_file().c_str());
  if (!ifs.fail() || ifs.good())
  {
    //throw InitFailedException("Cannot read spectrum "
      //  "from file " + datab.get_data_file());

    size_t i = 0;
    const size_t buf_len = 256;
    char buf[buf_len];
    while (ifs.good())
    {
      if (i == _eye_wl.size())
      {
        size_t n_new = _eye_wl.size() + 100;
        _eye_wl.reserve(n_new);
        _eye_value.reserve(n_new);
      }
      ifs.getline(buf, buf_len);
      if (buf[0] != '#')
      {
        istringstream in(buf);
        double l, s;
        if (in >> l >> s)
        {
          _eye_value.push_back(s);
          // conversion from nm^-1 to J
          _eye_wl.push_back(l);
          i++;
        }
      }
    }
    ifs.close();
  }
  else
    Messages::warning("Couldn't find (AVT spectrum)");


  Database db_spec;
  ifstream ifs_spect;
  db_spec.set_material("Emission-spectrum", get_option("Emission-spectrum-file", ""));
  ifs_spect.open(db_spec.get_data_file().c_str());
  if (!ifs_spect.fail() || ifs_spect.good())
  {
    //throw InitFailedException("Cannot read spectrum "
      //  "from file " + db_spec.get_data_file());

    size_t i = 0;
    const size_t buf_len = 256;
    char buf[buf_len];
    while (ifs_spect.good())
    {
      if (i == _eye_wl.size())
      {
        size_t n_new = _eye_wl.size() + 100;
        _emission_wl.reserve(n_new);
        _emission_value.reserve(n_new);
      }
      ifs_spect.getline(buf, buf_len);
      if (buf[0] != '#')
      {
        istringstream in(buf);
        double l, s;
        if (in >> l >> s)
        {
          _emission_value.push_back(s);
          // conversion from nm^-1 to J
          _emission_wl.push_back(l);
          i++;
        }
      }
    }
    ifs.close();
  }
  else
    Messages::warning("Couldn't find Internal source spectrum");
}


void
Tmm::do_setup_solution_variables(void)
{
  // we declare our solution variables
  declare_solution(GenerationRate, REAL, CELL, "1/cm^3/s");
  declare_solution(Intensity, REAL, CELL, "W/m^2");
  declare_solution(External_Source_ElectricField, REAL, CELL, "V/m");
  declare_solution(External_Source_Poynting, REAL, CELL, "W/m^2");

  declare_solution(Internal_Source_ElectricField, REAL, CELL, "a.u");
  declare_solution(Internal_Poynting, REAL, CELL, "a.u.");
  declare_solution(Internal_Power, REAL, CELL, "a.u.");
  declare_solution(Internal_Absorption, REAL, CELL, "a.u.");
  declare_solution(Internal_Intensity, REAL, CELL, "W/m^2");

  declare_solution(Energy_Loss, REAL, CELL, "W/m^3");


  declare_solution(Transmission, REAL, CELL, "1");
  declare_solution(Reflection, REAL, CELL, "1");
  declare_solution(Absorption, REAL, CELL, "1");

  declare_solution(AVT, REAL, CELL, "%");
  declare_solution(Generation_regions, REAL, CELL, "%");
  declare_solution(Polar, REAL, CELL, "a.u");


}


void
Tmm::do_solve(void)
{
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  NumericVector<libMesh::Number>& solution = system.get_local_solution_vector();
  solution.close();
  solution.zero();
  double dipole_sim_done = 0;
  reset_global_variables();
  string outdir = get_output_directory();
  //std::ifstream file_test(outdir + "/green_matrix.dat");
  //if (file_test.is_open())
  //  _green_vector_solved = 1;
  _green_vector_solved = 0;






  const MeshBase& mesh = get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  // NB: Tibercad by default uses length-scale in meters 
  //     This means that FEM derivatives d/dx are in 1/m
  //     To change this behavior it is necessary to define a different 'scaling'
  //     For instance if we want to use mesh_units in the assembly we need to:
  //     1. set the scaling to mesh units:
  //        get_scaling().set_length_scaling(get_mesh_units());
  //     2. use  build_finite_element(dim, fe_type, true)  
  //                                                ^ false is the default  
  //                                                
  // Now 2nd derivatives will be 1/mesh_units^2
  // We need a factor to transform rho/eps0 into V/mesh_units^2
  // Charge density is cm^-3, and Constants::e is in Coulomb, 
  // Constant::e0 is in C/Vm
  // The factor Lambda is such that rho*Lambda is in V/mesh_units^2
  // BUT (BUT) 
  // This is not that clever! Since Displacement and Polarization are already in C/m^2
  // it is easier to work with the derivatives in 1/m and rho/eps0 in V/m^2 
  // The factor 1e6 is for cm^3 -> m^3 in rho 
  get_scaling().set_length_scaling(1.0);

  DofMap& dof_map = system.get_dof_map();
  double mesh_units = get_mesh_units();
  vector<unsigned int> dof_indices;
  vector<double> intensity_integral;
  vector<double> Electric_Field_integral;
  vector<double> generation_rate_integral;
  vector<double> poynting_external_integral;
  vector<double> energy_loss_integral;
  vector<double> gen;
  vector<double> lambda_interp;
  vector<double> sun_interp;
  vector<double> Regions;

  vector<double> eye_num;
  vector<double> eye_dum;
  vector<double> eye_interp;

  vector<std::string> Names;

  vector<double> dipoles_power_global;

  double external_source_simulation{ 0 };
  double internal_source_simulation{ 0 };

  vector<complex<double>> Ki;
  if (_wavelength_vector.empty())
    for (size_t i = _down_lambda; i <= _up_lambda; i += _wavelength_steps)
      lambda_interp.push_back(i);
  else
    for (size_t i = 0; i < _wavelength_vector.size(); ++i)
      lambda_interp.push_back(_wavelength_vector[i]);
  if (!_spectrum.empty())
    sun_interp = Tmm::linear_interpolation1(_lambda, _spectrum, lambda_interp);
  if (!_eye_value.empty())
    eye_interp = Tmm::linear_interpolation1(_eye_wl, _eye_value, lambda_interp);

  Utils::Progress progress("Sweeping wavelength: ", lambda_interp.size());
  for (unsigned int i = 0; i < lambda_interp.size(); ++i)   //loop over wavelength
  {
    progress.progress_message();
    double eye{ 0 };
    if (!_eye_value.empty())
      eye = eye_interp[i];

    const double lambda = lambda_interp[i];
    _Wavelength.push_back(lambda);

    const unsigned int uvar = system.variable_number("G");

    // TODO this will not work if the 1D mesh is distributed. In that case, MPI calls could be used
    // to gather pieces from all processes
    MeshBase::const_element_iterator       el = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

    const double E = 1.0;
    const double incoming_angle = _incident_angle;
    const double c0 = 2.998e8 * 1e9;
    const double e0 = 8.85e-12;
    const double w = 2 * M_PI * c0 / lambda;
    const double plank_const = 1.055e-34;



    external_source_simulation = 0;
    internal_source_simulation = 0;


    vector<double> n_real;
    vector<double> n_imag;
    vector<double> l_length;
    vector<double> l;
    vector<double> Incoh;
    vector<double> elem_coordinate;;
    vector<double> eps_real;
    vector<double> eps_imag;
    vector<double> dipole_power;


    vector<int> BC_check;
    vector<int> BC_check_init;
    matrix2by2 Boundry_condition(1, 0, 0, 1);

    std::string direction;


    //**********************************************************************************************

    for (; el != end_el; ++el)
    {

      const Elem* elem = *el;

      dof_map.dof_indices(elem, dof_indices, uvar);
      const unsigned int n_dofs = dof_indices.size();
      TmmBulkModel& mod = *get_bulk_model<TmmBulkModel>(elem);
      mod.reinit(elem);


      l_length.push_back(dof_indices[0]);

      libMesh::Point pp = elem->vertex_average();
      mod.calculate(elem, pp(0), lambda);
      dipole_power.push_back(mod.get_emission_power());
      elem_coordinate.push_back(pp(0));

      if (i == 0)
      {
        ID sub_id = elem->subdomain_id();
        std::string name = get_environment().get_device().get_region_name(sub_id);
        Names.push_back(name);
      }


      libMesh::Complex nk = mod.get_refractive_index(lambda);
      Incoh.push_back(mod.get_coherent_index());
      n_real.push_back(real(nk));
      n_imag.push_back(imag(nk));
      eps_imag.push_back(2 * real(nk) * imag(nk)); //imaginary part of the epsilon
      eps_real.push_back(real(nk) * real(nk) - imag(nk) * imag(nk)); //imaginary part of the epsilon
      l.push_back(elem->volume() * (mesh_units / 1e-9));
      for (unsigned int s = 0; s < elem->n_sides(); s++)
      {
        TmmBoundaryModel* mod_int =
          get_interface_model<TmmBoundaryModel>(elem, s);
        if (mod_int != NULL)
        {
          if (mod_int->read_type() == "Mirror") {
            mod_int->Calculate_M_Matrix();
            BC_check_init.push_back(1);
            Boundry_condition.set(0, mod_int->get_element(0));
            Boundry_condition.set(1, mod_int->get_element(1));
            Boundry_condition.set(2, mod_int->get_element(2));
            Boundry_condition.set(3, mod_int->get_element(3));
          }
          if (mod_int->read_type() == "Incident Wave") {
            external_source_simulation = 1;
            BC_check_init.push_back(0);
            if (s % 2 == 0)
              direction = "top to down propagation";
            else
              direction = "down to top propagation";

          }
        }
        if (mod_int == NULL)
          BC_check_init.push_back(0);
      }
    }

    size_t mesh_size = n_real.size();

    if (direction == "top to down propagation") // if the incomming wave is from the other sild, so the vectors should be revesed
    {
      for (size_t uu = 1; uu < BC_check_init.size();uu += 2)
        BC_check.push_back(BC_check_init[uu]);
    }
    if (i == 0)
    {
      for (size_t lm = 0; lm < Names.size();lm++)
        if (lm == 0 || Names[lm] != Names[lm - 1])
        {
          _regions_name.push_back(Names[lm]);
          Regions.push_back(lm);
        }
      Regions.push_back(mesh_size);
      _Generation_rate.resize(mesh_size);
      _Poynting_external.resize(mesh_size);
      _Intensity.resize(mesh_size);
      _External_Source_ElectricField.resize(mesh_size);
      _Energy_loss_external.resize(mesh_size);

      _Internal_Source_ElectricField.resize(mesh_size);
      _Internal_Poynting.resize(mesh_size);
      _Internal_Power.resize(mesh_size);
      _Internal_Absorption.resize(mesh_size);
      _Abs.resize(mesh_size);
      _Internal_Intensity.resize(mesh_size);
      if (!_green_vector_solved)
        _green_vector.resize(mesh_size);  //resize if internal source hasn't been solved yet	 
    }

    //****************************************************************************
    //*****************************External_source_simulation******************************
    if (external_source_simulation)
    {

      if (sun_interp.empty())
        throw InitFailedException("External illumination spectrum"
          "is missing ");
      double radiation = sun_interp[i] / 1e3;            // [W/m^2/nm]
      double Esun2 = 2 * radiation / (c0 * 1e-9 * e0);   // [V^2/m^2/nm]
      vector<double> theta(mesh_size, 0);
      //theta = Tmm::theta_cal(n_real,incoming_angle);
      double rnd = 1;
      double phase_step;
      vector<double> Generation_rate_avg(mesh_size, 0);
      vector<double> poynting_external_avg(mesh_size, 0);
      vector<double> Intensity_avg(mesh_size, 0);
      vector<double> Electric_Field_avg(mesh_size, 0);
      vector<double> energy_loss_avg(mesh_size, 0);

      double avg_reflection = 0;
      double avg_transmission = 0;
      double avg_Absorption = 0;


      double coh_mod = _coh_mod;

      if (coh_mod)
      {
        for (size_t nm = 0; nm < Incoh.size(); ++nm)
          if (Incoh[nm] == 1)
            rnd = _coh_mod;
      }
      else
      {
        for (size_t nm = 0; nm < Incoh.size(); ++nm)
          if (Incoh[nm] == 1)
            rnd = 5;
      }

      phase_step = 2 * M_PI / rnd;
      double N0 = 0;
      //***********************************************************************
        //**********loop added for highly absorbing devices*******************
      // check for higly absorbing material
        // if this is the case, we will iqnore N0 number elements	from the back.  
      matrix2by2 Unit(1, 0, 0, 1);
      matrix2by2 DD(1, 0, 0, 1);
      matrix2by2 MM(0, 0, 0, 0);
      matrix2by2 TT_load(0, 0, 0, 0);
      matrix2by2 TT(1, 0, 0, 1);
      for (size_t k = 0;k <= mesh_size - 2; ++k)
      {
        DD = get_D(n_real[k], n_imag[k], n_real[k + 1], n_imag[k + 1], 0, lambda, 0);
        if (Incoh[k] == 1 && Incoh[k + 1] == 0)
        {
          MM = get_M(n_real[k], n_imag[k], l[k], lambda, 0, 0);
        }
        else
        {
          MM = get_M(n_real[k], n_imag[k], l[k], lambda, 0, 0);
        }
        TT_load = MM * DD;
        if (real(abs(TT.get(0))) > 1e150)
        {
          TT_load = Unit;
          TT = TT * TT_load;
          if (N0 == 0)
            N0 = mesh_size - k;
        }
        else
          TT = TT * TT_load;
      }

      srand(time(NULL));
      for (size_t iter = 0; iter < rnd; ++iter)    // loop over added phases
      {
        double phase;
        double r = _reflectivity;

        //***********************************************************************
        //**********************Main TMM loop************************************
        matrix2by2 T(1, 0, 0, 1);
        matrix2by2 E_N(1, 0, -r, 0);
        matrix2by2 E_I(0, 0, 0, 0);
        vector<complex<double>> E_F;
        vector<complex<double>> E_B;
        if (N0 != 0)
          for (size_t nm = 0; nm < N0; nm++)
          {
            E_F.push_back(0);
            E_B.push_back(0);
          }
        E_F.push_back(1);
        E_B.push_back(-r);
        //T = get_M(n_real[n_real.size()-N0-2],n_imag[n_real.size()-N0-2],l[n_real.size()-N0-2],lambda,theta[n_real.size()-N0-2],0);
        for (int k = mesh_size - N0 - 2;k >= 0; --k)
        {
          DD = get_D(n_real[k], -n_imag[k], n_real[k + 1], -n_imag[k + 1], 0, lambda, 0);
          if (Incoh[k] == 1 && Incoh[k + 1] == 0)
          {
            if (coh_mod)
            {
              phase = 2 * M_PI * (double)rand() / RAND_MAX;
            }
            else
            {
              phase = phase_step * iter;
            }
            MM = get_M(n_real[k], n_imag[k], l[k], lambda, 0, phase);
          }
          else
          {
            MM = get_M(n_real[k], n_imag[k], l[k], lambda, 0, 0);
          }
          TT_load = DD * T;
          T = MM * TT_load;
          E_I = T * E_N;
          E_F.push_back(E_I.get(0));
          E_B.push_back(E_I.get(2));
        }
        //****************************************************************************
        //***************normalizing electric field matrix****************************
        vector<complex<double>> E_F_NORM(mesh_size);
        vector<complex<double>> E_B_NORM(mesh_size);
        for (size_t nm = 0; nm < mesh_size; ++nm)
        {
          E_F_NORM[nm] = E_F[E_F.size() - 1 - nm] / E_F[E_F.size() - 1];
          E_B_NORM[nm] = E_B[E_F.size() - 1 - nm] / E_F[E_F.size() - 1];
        }
        //***********************************************************************
        //**********reflection and transmission calculation**********************
        complex<double> Refl, Trans;
        Refl = pow(abs(T.get(2) / T.get(0)), 2);
        avg_reflection = avg_reflection + real(Refl) / rnd;
        complex<double> nc_first(n_real[0], n_imag[0]);
        complex<double> nc_last(n_real[mesh_size - 1], n_imag[n_imag.size() - 1]);
        complex<double> ratio_complex;
        ratio_complex = ((nc_last)*cos(theta[theta.size() - 1] * M_PI / 180)) / (nc_first * cos(theta[0] * M_PI / 180));
        Trans = ratio_complex * pow(abs(1.0 / T.get(0)), 2);
        avg_transmission = avg_transmission + real(Trans) / rnd;
        avg_Absorption = avg_Absorption + (1 - real(Refl) - real(Trans)) / rnd;

        //****************************************************************************
        //***************Calculating average values over random phases****************
        complex<double> Etot;
        complex<double> H_F;
        complex<double> H_B;
        complex<double> Intensity;
        double Generation_rate;
        double poynting_external;
        double energy_loss;
        for (size_t nm = 0; nm < mesh_size; ++nm)
        {
          Etot = E_F_NORM[nm] + E_B_NORM[nm]; //[V/m]
          Electric_Field_avg[nm] += (real(Etot) * sqrt(Esun2) / rnd); // [V/m/nm^0.5] 
          complex<double> N_i(n_real[nm], n_imag[nm]);

          H_F = 1e-9 * c0 * e0 * lambda * N_i / (2 * M_PI) * E_F_NORM[nm];
          H_B = 1e-9 * c0 * e0 * lambda * N_i / (2 * M_PI) * E_B_NORM[nm];

          poynting_external = (0.5 * real((E_F_NORM[nm] * conj(H_F)) - (E_B_NORM[nm] * conj(H_B))));
          poynting_external_avg[nm] = poynting_external_avg[nm] + poynting_external / rnd; //[W/m^2/nm]

          Intensity = (0.5 * c0 * 1e-9 * e0 * n_real[nm] * Esun2 * pow(abs(Etot), 2)); //[W/m^2/nm]
          Intensity_avg[nm] = Intensity_avg[nm] + real(Intensity) / rnd;
          Generation_rate = real(1 / (plank_const * w) * (4 * M_PI * n_imag[nm] * 1e7 / (lambda)) * real(Intensity) / 1e4); //[1/cm^3/s/nm]
          Generation_rate_avg[nm] += Generation_rate / rnd;
          energy_loss = (M_PI * c0 / lambda * eps_imag[nm] * Esun2 * pow(abs(Etot), 2));
          energy_loss_avg[nm] = energy_loss_avg[nm] + energy_loss / rnd;
        }
      }
      vector<double> sumation;
      for (size_t re = 1; re < Regions.size(); re++)
      {
        double load = 0;
        for (size_t el = Regions[re - 1]; el < Regions[re]; el++)
        {
          load += Generation_rate_avg[el];
        }
        sumation.push_back(load);
      }
      double load = 0;
      for (size_t re = 0; re < sumation.size();re++)
        load += sumation[re];
      for (size_t re = 0; re < sumation.size();re++)
        sumation[re] = (sumation[re] / load) * avg_Absorption;

      _Generation_regions.push_back(sumation);
      _Electric_Field_External.push_back(Electric_Field_avg);   // [V/m]

      const double scale = 1e-6;
      _Reflection.push_back((int)(real(avg_reflection) / scale) * scale);
      _Transmission.push_back((int)(real(avg_transmission) / scale) * scale);
      _Absorption.push_back((int)(avg_Absorption / scale) * scale);

      //****************************************************************************
      //*****************************AVT Calculation********************************
      if (!_eye_value.empty())
      {
        eye_num.push_back(eye * radiation * _Transmission[_Transmission.size() - 1]);
        eye_dum.push_back(eye * radiation);
      }
      //****************************************************************************
      //*******************Calculating integral over wavelengths********************
      double delta_wl;
      double delta_f;
      if (i == 0) // first wavelength
      {
        if (lambda_interp.size() > 1)
        {
          delta_wl = abs(lambda_interp[i + 1] - lambda_interp[i]) / 2; //first wavlength
          delta_f = abs((c0 / lambda_interp[i + 1]) - (c0 / lambda_interp[i])) / 2;
        }
        else
        {
          delta_wl = 1;
          delta_f = 1;
        }

        Electric_Field_integral.resize(mesh_size);
        intensity_integral.resize(mesh_size);
        generation_rate_integral.resize(mesh_size);
        poynting_external_integral.resize(mesh_size);
        energy_loss_integral.resize(mesh_size);
      }
      else
      {
        if (i == lambda_interp.size() - 1)
        {
          delta_wl = abs(lambda_interp[i] - lambda_interp[i - 1]) / 2; //last wavelength    // [nm]
          delta_f = abs((c0 / lambda_interp[i]) - (c0 / lambda_interp[i - 1])) / 2;
        }
        else
        {
          delta_wl = (lambda_interp[i + 1] - lambda_interp[i - 1]) / 2; //midle wavelengths
          delta_f = abs((c0 / lambda_interp[i + 1]) - (c0 / lambda_interp[i - 1])) / 2;
        }
      }



      for (size_t km = 0; km < mesh_size; km++)
      {
        Electric_Field_integral[km] += Electric_Field_avg[km] * sqrt(delta_wl);  // [V/m]
        intensity_integral[km] += Intensity_avg[km] * delta_wl;   // [W/m^2]
        generation_rate_integral[km] += Generation_rate_avg[km] * delta_wl;   // [1/cm^3/s]
        poynting_external_integral[km] += poynting_external_avg[km] * delta_wl;   // [W/m^2]
        energy_loss_integral[km] += energy_loss_avg[km] * delta_wl;
      }


      if (i == lambda_interp.size() - 1)                   //last wavelength
        for (size_t nm = 0; nm < l_length.size(); ++nm)
        {
          _Generation_rate[nm] = generation_rate_integral[nm];
          _Poynting_external[nm] = poynting_external_integral[nm];
          _Intensity[nm] = intensity_integral[nm];
          _External_Source_ElectricField[nm] = Electric_Field_integral[nm];
          _Energy_loss_external[nm] = energy_loss_integral[nm];
          if (!_eye_value.empty())
          {
            _AVT = 0;
            double num = 0;
            double dum = 0;
            for (size_t mm = 0; mm < eye_num.size();mm++)
            {
              num += eye_num[mm];
              dum += eye_dum[mm];    ///asumming 1nm sweep of wavelength
            }
            _AVT = num / dum * 100;
          }

        }


    }
    //*****************************End of External_source_simulation******************************

    //****************************************************************************
    //*****************************Internal_source_simulation******************************
    vector<double> dipole_loc_list;
    vector<double> dipole_pow_list;
    if (!_dipole_coordinate.empty()) // if dipole deffined as list of points, not defined by dd module
      for (double dipole_num = 0;dipole_num < _dipole_coordinate.size();++dipole_num)
        for (double elem = 1; elem < elem_coordinate.size(); ++elem)
          if (_dipole_coordinate[dipole_num] >= elem_coordinate[elem - 1] && _dipole_coordinate[dipole_num] < elem_coordinate[elem])
          {
            dipole_loc_list.push_back(elem - 1); // elements contains dipole
            dipole_pow_list.push_back(_dipole_power[dipole_num]); // power of each dipole
            internal_source_simulation = 1; // solve the internal source
          }

    double dipole_wl_index; //index of lambda in emiison speectrum file
    if (!dipole_power.empty() && !_emission_value.empty())  // if internal source define based on coupling with dd module
    {
      for (size_t elem = 0; elem < dipole_power.size(); ++elem)
        if (abs(dipole_power[elem]) >= 1e-18)      // minimum intensity to solve for 
        {
          dipole_loc_list.push_back(elem); // elements contains dipole
          dipole_pow_list.push_back(dipole_power[elem]); // power of each dipole
          for (size_t el = 0;el < _emission_wl.size();el++)
            if (lambda == int(_emission_wl[el]) && _emission_value[el] >= 0.001)
            {
              internal_source_simulation = 1;  // solve the internal source
              dipole_wl_index = el;
            }
        }

    }
    //store the emission rate in vector to avoid repeating simulation
    // of teh internal source.
    if (dipoles_power_global.empty() && internal_source_simulation)
      for (size_t nm = 0; nm < dipole_pow_list.size(); nm++)
        dipoles_power_global.push_back(dipole_pow_list[nm]);



    if (internal_source_simulation && !_green_vector_solved)  //solve if there is intrnal emission and the _green_vector is not ready yet
    {
      //std::cout<<"start solving ... "<<std::endl;
      if (_emission_value.empty())
        throw InitFailedException("Internal illumination spectrum"
          "is missing ");

      _Internal_Wavelength.push_back(lambda);


      vector<double> kr_vec;
      std::vector<double> outcoupling;

      double ks = 2 * M_PI * n_real[dipole_loc_list[0]] / lambda;
      if (_steps != 0)
      {
        kr_vec.resize(_steps);
        kr_vec[0] = 0;
        for (size_t u = 1; u < _steps; ++u)
        {
          kr_vec[u] = kr_vec[u - 1] + ks * _ratio[0] / _steps;
        }
      }
      else
      {
        _steps = 1;
        kr_vec.resize(_steps);
        kr_vec[0] = ks * _ratio[0];
      }

      double A_P;
      double A_N;

      double rnd = 1;
      double phase_step = 0;

      _angle.resize(_steps);
      _Poynting_front.resize(_steps);
      _Poynting_back.resize(_steps);

      for (size_t dipole_num = 0;dipole_num < dipole_loc_list.size();++dipole_num)
      {
        double dipole_loc = dipole_loc_list[dipole_num];
        vector<double> internal_power_per_dipole(mesh_size, 0);
        for (auto Oraintation : _oraintation_vec)
        {


          //*****************************loop for number of dipoles******************************
          for (auto Mode : _polarization_vec)
          {

            std::cout << "Solving for internal source --> dipole element is " << dipole_loc << " power is : " << dipole_pow_list[dipole_num] << " lambda is  " << lambda << std::endl;
            complex<double> Ns(n_real[dipole_loc], n_imag[dipole_loc]);
            // std::cout<<"Ns is : "<<Ns <<std::endl;
                 //*****************************loop for radial wave numbers******************************
            double cnt = 0;
            double kr;
            for (size_t kr_num = 0; kr_num < _steps; ++kr_num)
            {
              kr = kr_vec[kr_num];
              if (dipole_num == 0)
              {
                for (size_t nm = 0; nm < dipole_power.size(); ++nm)
                  _green_vector[nm].resize(dipole_pow_list.size());
              }

              complex<double> ki_s((2 * M_PI * n_real[dipole_loc]) / lambda, (2 * M_PI * n_imag[dipole_loc]) / lambda);
              double cos_phi_inter;
              cos_phi_inter = real(cmlx_sqrt(complex<double>(1.0 - (pow(kr, 2) / pow(ki_s, 2)))));
              if (cos_phi_inter > 1)
                cos_phi_inter = 1;
              _angle[cnt] = acos(cos_phi_inter);

              // Calculating Source terms with respect to dipole properties
              dipole_source(A_P, A_N, Mode, Oraintation, lambda, cos_phi_inter);

              // std::cout<< dipole_pow_list[dipole_num]<< "		"<<A_P<<std::endl;
                  //*****************************calculating electric field of two ends of the device******************************
              vector<complex<double>> E_int_f(mesh_size);
              vector<complex<double>> E_int_b(mesh_size);
              vector<complex<double>> E_int(mesh_size);
              vector<complex<double>> E_int_norm(mesh_size);

              for (size_t iter = 0; iter < rnd; ++iter)    // loop over added phases
              {
                solving_internal_source(E_int_f, E_int_b, n_real, n_imag, l, lambda, dipole_loc, kr, Mode, A_P, phase_step * iter);
                //std::cout<<"phase_step * iter" << "	 =	"<< phase_step * iter<<std::endl;
                for (size_t nm = 0; nm < E_int_f.size(); ++nm)
                  E_int[nm] += (E_int_f[nm] + E_int_b[nm]) / rnd;
              }


              double coefficient = 1 / _steps / _oraintation_vec.size() / _polarization_vec.size();  ///lambda_interp.size()
              
              for (size_t nm = 0; nm < E_int_f.size(); ++nm)
              {
                _Internal_Source_ElectricField[nm] += real(E_int[nm]) * coefficient;
                _Internal_Intensity[nm] += 0.5 * c0 * 1e-9 * e0 * n_real[nm] * pow(abs(E_int[nm]), 2);
                // std::cout<<nm << "		"<< _Internal_Intensity[nm]<<std::endl;
              }


              if (kr_num == 0 && dipole_num == 0 && Oraintation == 0 && Mode == 0)
              {
                std::vector<double> buffer(mesh_size);
                for (size_t cnter = 0; cnter < mesh_size;cnter++)
                  buffer[cnter] = real(E_int[cnter]);
                _Electric_Field_Internal.push_back(buffer);
              }


              //**************a loop to calculate Poynting Vector from electric field******************************
              double poynting;
              double internal_power;
              complex<double> Ns(n_real[dipole_loc], n_imag[dipole_loc]);
              complex<double> kis((2 * M_PI * n_real[dipole_loc]) / lambda, (2 * M_PI * n_imag[dipole_loc]) / lambda);
              complex<double> kzs;
              kzs = kis * cmlx_sqrt(complex<double>(1.0 - (pow(kr, 2) / pow(kis, 2))));
              complex<double> kzs2;
              kzs2 = pow(kzs, 2);
              double dkr = ks * _ratio[0] / _steps;
              for (size_t nm = 0; nm < mesh_size;nm++)
              {
                complex<double> ki((2 * M_PI * n_real[nm]) / lambda, (2 * M_PI * n_imag[nm]) / lambda);
                complex<double> kzi(1, 0);

                kzi = ki * cmlx_sqrt(complex<double>(1.0 - (pow(kr, 2) / pow(ki, 2))));
                complex<double> cos_phi;
                cos_phi = kzi / ki;
                complex<double> Ni(n_real[nm], n_imag[nm]);
                if (Mode == 0)  //TE Mode
                  poynting = real(Ni * cos_phi * conj(E_int_f[nm] + E_int_b[nm]) * (E_int_f[nm] - E_int_b[nm]));
                else    //TM Mode
                  poynting = real(Ni * conj(cos_phi) * conj(E_int_f[nm] + E_int_b[nm]) * (E_int_f[nm] - E_int_b[nm]));

                _Internal_Poynting[nm] += poynting * coefficient;
                internal_power = real((2 * M_PI) * (poynting * coefficient) * kr * dkr / Ns / kzs2);
                _Internal_Power[nm] += internal_power;
                internal_power_per_dipole[nm] += poynting * coefficient;

              }

              //*****************************Assigning Output Power of two ends of the device ******************************

              vector<double> Abs(mesh_size);
              double absorb = 1;

              double sum = 0;
              for (size_t nm = 1; nm < mesh_size; ++nm)
              {
                Abs[nm] = 0.5 * c0 * 1e-9 * e0 * n_real[nm] * pow(abs(E_int[nm]), 2) * n_imag[nm];
                sum += Abs[nm];
              }

              for (size_t nm = 1; nm < mesh_size; ++nm)
              {
                if (sum != 0)
                  Abs[nm] = Abs[nm] * absorb * coefficient * _emission_value[dipole_wl_index] / sum;
              }

              for (size_t nm = 0; nm < Abs.size(); ++nm)
              {
                _Internal_Absorption[nm] += Abs[nm] * dipole_pow_list[dipole_num];
                _green_vector[nm][dipole_num] += Abs[nm];
                dipole_sim_done = 1;
              }


              ++cnt;

            } // end of radial wave numbers loop

          } // end of number of modes loop

        } // end of number of orain loop
        double P0, Pend, Pdip, extraction;
        P0 = abs(internal_power_per_dipole[0]);
        Pend = abs(internal_power_per_dipole[mesh_size - 1]);


        Pdip = abs(internal_power_per_dipole[dipole_loc]) + abs(internal_power_per_dipole[dipole_loc - 2]);

        if (Pdip != 0)
        {
          extraction = (P0 + Pend) / Pdip;
        }
        else
        {
          extraction = 0;
        }
        double sum_power = 0;
        for (auto dipole_power : dipole_pow_list)
          sum_power += dipole_power;
        outcoupling.push_back(extraction * dipole_pow_list[dipole_num]);
      }// end of number of dipoles
      double avg = 0;

      for (auto dipole_outcoupling : outcoupling)
        avg += dipole_outcoupling;


      _OutCoupling.push_back(avg / outcoupling.size()); //pushing the average outcoupling for dipoles


      vector<double> sumation;
      for (size_t re = 1; re < Regions.size(); re++)
      {
        double load = 0;
        for (size_t el = Regions[re - 1]; el < Regions[re]; el++)
          load += _Internal_Absorption[el];
        sumation.push_back(load);
      }
      double load = 0;
      for (size_t re = 0; re < sumation.size();re++)
        load += sumation[re];

      for (size_t re = 0; re < sumation.size();re++)
        if (load != 0)
          sumation[re] = (sumation[re] / load);
        else
          sumation[re] = 0;
      _Generation_regions_internal.push_back(sumation);


    }// end of dipole simulation




  }// end of loop of wavelength


  std::vector<double> Internal_Absorption_green(_Internal_Absorption.size());

  if (dipole_sim_done && !_green_vector_solved) // store _green_vector in file
  {
    std::cout << " \033[33m saving green_matrix data  \033[0m " << std::endl;
    _green_vector_solved = 1;
    std::ofstream outFile(outdir + "/green_matrix.dat");
    if (!outFile)
      std::cout << "\033[33m Cannot creat green_matrix.dat file!!! \033[0m " << std::endl;


    size_t rows = _green_vector.size();
    size_t cols = _green_vector[0].size();
    outFile << rows << " " << cols << std::endl;

    // Now write the _green_vector elements
    for (const auto& row : _green_vector) {
      for (size_t j = 0; j < row.size(); ++j) {
        outFile << row[j]; // Write the value as it is
        if (j < row.size() - 1) {
          outFile << " "; // Separate values by space
        }
      }
      outFile << std::endl; // Newline after each row
    }

    outFile.close();


  }
  else if (_green_vector_solved && !dipoles_power_global.empty())  // read from file to complet the _green_vector in the next iteration
  {

    //std::cout<<"Optical simulation is over-2"<<std::endl;  
    std::cout << " \033[33m reading and completting green_matrix data  \033[0m " << std::endl;
    _green_vector.clear();
    std::ifstream inFile(outdir + "/green_matrix.dat");
    if (!inFile)
    {
      std::cout << "\033[33m Cannot open green_matrix.dat file! \033[0m " << std::endl;
    }
    size_t rows, cols;
    inFile >> rows >> cols;
    _green_vector.resize(rows);
    //std::cout<<"Optical simulation is over-1"<<std::endl;
      // Create the _green_vector
    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < cols; ++j) {
        _green_vector[i].resize(cols);
        inFile >> _green_vector[i][j]; // Read each element
      }
    }
    inFile.close();

  }
  //std::cout<<"Optical simulation is over0"<<std::endl;

  if (!dipoles_power_global.empty() && _green_vector_solved)
    for (size_t i = 0; i < _green_vector.size(); ++i)
      for (size_t j = 0; j < _green_vector[i].size(); ++j)
      {
        Internal_Absorption_green[i] += _green_vector[i][j] * dipoles_power_global[j];
        _Internal_Absorption[i] = Internal_Absorption_green[i];
      }

  for (size_t nm = 0; nm < _Generation_rate.size(); ++nm)  // suming generation rate for internal and external sources
    _Generation_rate[nm] += Internal_Absorption_green[nm];
  std::cout << "Optical simulation is over" << std::endl;


  solution.close();
  system.update();
}



void
Tmm::do_print_info(void)
{
  Messages::info("1D TMM");
}


PhysicalModel*
Tmm::create_bulk_model(const ModelOptions& options,
  const Material* mat) const
{
  return TmmBulkModel::create(mat, options);
}



PhysicalModel*
Tmm::create_boundary_model(const ModelOptions& options,
  const MaterialBoundary* boundary) const
{
  return(TmmBoundaryModel::create(boundary, options));
}

void
Tmm::plot_globaldata(void)
{
  //std::cout<<"Tmm::plot_globaldata(void) "<<std::endl;
  string outdir = get_output_directory();
  //std::cout<<"AVT IS "<<_solutions.count(AVT)<<std::endl;
  if (_solutions.count(AVT))
    if (_AVT != 0)
    {
      string filename(outdir + "/" + get_output_filename() + "_AVT.dat");
      ofstream file;
      file.open(filename.c_str());
      if (file.good())
      {
        file << "# " << get_type() << " AVT (" << get_name() << ")\n";
        file << "# " << 1 << " AVT" << "\n";
        file << "# " << "AVT" << "\n";

        file << _AVT << " ";
        file << "\n";

      }


      file.close();

    }

  if (_solutions.count(Generation_regions))
    if (!_Generation_regions.empty())
    {


      string filename(outdir + "/" + get_output_filename() + "_Absorption_Regions.dat");

      ofstream file;
      file.open(filename.c_str());

      if (file.good())
      {
        file << "# " << get_type() << " Absorption_Regions (" << get_name() << ")\n";
        file << "# " << 0 << " WaveLength [nm]" << "\n";
        for (unsigned int i = 0; i < _regions_name.size(); ++i)
        {
          file << "# " << i + 1 << " Region => " << _regions_name[i] << "\n";
        }

        file << "# " << "WaveLength ";
        for (unsigned int i = 0; i < _regions_name.size(); ++i)
        {
          file << _regions_name[i] << " ";
        }
        file << "\n";

        for (unsigned int i = 0; i < _Wavelength.size(); i++)
        {
          file << _Wavelength[i] << " ";
          for (unsigned int j = 0; j < _regions_name.size(); j++)
            file << _Generation_regions[i][j] << " ";
          file << "\n";
        }
      }


      file.close();
    }

  if (_solutions.count(Generation_regions))

    if (!_Generation_regions_internal.empty())
    {

      //std::cout<<"if (!_Generation_regions_internal.empty()) "<<std::endl;


      string filename(outdir + "/" + get_output_filename() + "_Absorption_Regions_Internal.dat");

      ofstream file;
      file.open(filename.c_str());

      if (file.good())
      {
        file << "# " << get_type() << " Absorption_Regions_Internal (" << get_name() << ")\n";
        file << "# " << 0 << " Wavelength" << "\n";
        file << "# " << 1 << " OutCoupling" << "\n";
        for (unsigned int i = 0; i < _regions_name.size(); ++i)
        {
          file << "# " << i + 2 << " Region => " << _regions_name[i] << "\n";
        }

        file << "# " << "Wavelength " << "OutCoupling ";
        for (unsigned int i = 0; i < _regions_name.size(); ++i)
        {
          file << _regions_name[i] << " ";
        }
        file << "\n";

        for (unsigned int i = 0; i < _Generation_regions_internal.size(); i++) //internal emission Wavelength
        {
          file << _Internal_Wavelength[i] << " ";
          file << _OutCoupling[i] << " ";
          for (unsigned int j = 0; j < _regions_name.size(); j++)
            file << _Generation_regions_internal[i][j] << " ";
          file << "\n";
        }
      }


      file.close();
    }
  //std::cout<<"f (_solutions.count(External_Source_ElectricField))"<<std::endl;
  if (_solutions.count(External_Source_ElectricField))
    if (!_Electric_Field_External.empty())
    {
      string filename(outdir + "/" + get_output_filename() + "_electric_field_external.dat");

      ofstream file;
      file.open(filename.c_str());

      if (file.good())
      {
        unsigned int step = 1;
        std::vector<double> wl_vect;
        if (_Wavelength.size() < 10)
        {
          wl_vect.resize(_Wavelength.size());
          for (unsigned int i = 0; i < _Wavelength.size(); i++)
            wl_vect[i] = _Wavelength[i];
        }
        else
        {
          wl_vect.resize(10);
          step = (_Wavelength.size()) / 9;
          for (unsigned int i = 0; i < 10; i++)
            wl_vect[i] = _Wavelength[step * i];

        }

        // header
        file << "# " << get_type() << " TMMM (" << get_name() << ")\n";
        file << "# " << "Electric Field [V/m/nm]" << "\n";
        for (unsigned int i = 0; i < wl_vect.size(); ++i)
        {
          file << "# " << i << " WaveLength = " << wl_vect[i] << "[nm]" << "\n";
        }
        file << "# " << "x ";
        for (unsigned int i = 0; i < wl_vect.size(); ++i)
        {
          file << wl_vect[i] << "[nm] ";
        }
        file << "\n";

        for (unsigned int j = 0; j < _External_Source_ElectricField.size(); ++j)
        {
          file << j << " ";
          for (unsigned int i = 0; i < wl_vect.size(); ++i)
          {
            file << _Electric_Field_External[i * step][j] << " ";
          }
          file << "\n";
        }
      }
      file.close();
    }
  //std::cout<<" if (_solutions.count(Internal_Source_ElectricField))"<<std::endl;
  if (_solutions.count(Internal_Source_ElectricField))
    if (!_Electric_Field_Internal.empty())
    {
      string filename(outdir + "/" + get_output_filename() + "_electric_field_inernal.dat");

      ofstream file;
      file.open(filename.c_str());

      if (file.good())
      {
        unsigned int step = 1;
        std::vector<double> wl_vect;
        if (_Internal_Wavelength.size() < 10)
        {
          wl_vect.resize(_Internal_Wavelength.size());
          for (unsigned int i = 0; i < _Internal_Wavelength.size(); i++)
            wl_vect[i] = _Internal_Wavelength[i];
        }
        else
        {
          wl_vect.resize(10);
          //std::cout<<_Internal_Wavelength.size() << " vs " << floor(_Internal_Wavelength.size() / 9.0) <<std::endl;
          step = (_Internal_Wavelength.size()) / 9;
          for (unsigned int i = 0; i < 10; i++)
          {
            wl_vect[i] = _Internal_Wavelength[step * i];
            //std::cout<<i << " vs " <<step * i <<std::endl;
          }

        }
        
        file << "# " << get_type() << " TMM (" << get_name() << ")\n";
        for (unsigned int i = 0; i < wl_vect.size(); ++i)
        {
          file << "# " << i << " WaveLength = " << wl_vect[i] << "[nm]" << "\n";
        }
        file << "# " << "x ";
        for (unsigned int i = 0; i < wl_vect.size(); ++i)
        {
          file << wl_vect[i] << "[nm] ";
        }
        file << "\n";
        
        for (unsigned int j = 0; j < _Electric_Field_Internal.size(); ++j)
        {
          file << j << " ";
          for (unsigned int i = 0; i < wl_vect.size(); ++i)
          {
            file << _Electric_Field_Internal[j][i * step] << " ";
          }
          file << "\n";
        }
      }
      file.close();
    }

  //std::cout<<" if (_solutions.count(Transmission) || _solutions.count(Reflection) || _solutions.count(Absorption))"<<std::endl;

  if (_solutions.count(Transmission) || _solutions.count(Reflection) || _solutions.count(Absorption))
    if (!_Electric_Field_External.empty())
    {
      string filename(outdir + "/" + get_output_filename() + ".dat");

      ofstream file;
      file.open(filename.c_str());

      if (file.good())
      {
        //std::cout<<"start printing "<<std::endl;
          // header
        file << "# " << get_type() << " TMM (" << get_name() << ")\n";

        file << "# WaveLength[nm] " << " Transmission[1] "
          << " Reflection[1] "
          << " Absorption[1] " << "\n";

        for (unsigned int i = 0; i < _Transmission.size(); i++)
        {
          file << _Wavelength[i] << " "
            << _Transmission[i] << " "
            << _Reflection[i] << " "
            << _Absorption[i] << "\n";
        }
      }
      file.close();
    }

  //std::cout<<" END OF GLOBAL PLOTTING"<<std::endl;
}

void
Tmm::get_solution_secure(const Elem* elem,
  std::map<ID, std::vector<double> >& solutions,
  const std::vector<Point>& p)
{
  _solutions = solutions;
  unsigned int np = p.size();

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();
  const NumericVector<libMesh::Number>& solution = system.get_solution_vector();
  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system.get_dof_map();

  const unsigned int u_var = system.variable_number("G");

  FEType fe_type = system.variable_type(u_var);
  std::unique_ptr<FEBase> fe(build_finite_element(dim, fe_type));

  vector<unsigned int> dof_indices;

  // element shape functions
  const vector<vector<Real> >& phi = fe->get_phi();
  const vector<vector<RealGradient> >& dphi = fe->get_dphi();
  const vector<Point>& real_pts = fe->get_xyz();

  ID subdomain = elem->subdomain_id();

  fe->reinit(elem, &p);

  dof_map.dof_indices(elem, dof_indices, u_var);
  const unsigned int n_dofs = dof_indices.size();

  // cell data variable
  RealGradient field(0);
  TmmBulkModel& mod = *get_bulk_model<TmmBulkModel>(elem);

  if (solutions.count(Internal_Intensity))
  {
    solutions[Internal_Intensity][0] = _Internal_Intensity[dof_indices[0]];
  }

  if (solutions.count(Internal_Absorption))
  {
    solutions[Internal_Absorption][0] = _Internal_Absorption[dof_indices[0]];
  }

  if (solutions.count(Internal_Poynting))
  {
    solutions[Internal_Poynting][0] = _Internal_Poynting[dof_indices[0]];
  }

  if (solutions.count(Internal_Power))
  {
    solutions[Internal_Power][0] = _Internal_Power[dof_indices[0]];
  }

  if (solutions.count(Internal_Source_ElectricField))
  {
    solutions[Internal_Source_ElectricField][0] = _Internal_Source_ElectricField[dof_indices[0]];
  }

  if (solutions.count(External_Source_ElectricField))
  {
    solutions[External_Source_ElectricField][0] = _External_Source_ElectricField[dof_indices[0]];
  }

  if (solutions.count(Energy_Loss))
  {
    solutions[Energy_Loss][0] = _Energy_loss_external[dof_indices[0]];
  }
  if (solutions.count(Intensity))
  {
    solutions[Intensity][0] = _Intensity[dof_indices[0]];
  }

  if (solutions.count(GenerationRate))
  {
    solutions[GenerationRate][0] = _Generation_rate[dof_indices[0]];
  }
  if (solutions.count(External_Source_Poynting))
  {
    solutions[External_Source_Poynting][0] = _Poynting_external[dof_indices[0]];
  }
}


