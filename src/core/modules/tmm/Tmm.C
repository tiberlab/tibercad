// $Id$

#include "Tmm.h"
#include "TiberLinearSystem.h"
#include "TmmBulkModel.h"
#include "TmmBoundaryModel.h"
#include "SimulationEnvironment.h"
#include "Device.h"
#include "Messages.h"
#include "Database.h"
#include <cctype>


#include "libmesh/dof_map.h"




// This is needed in order to create the shared module library
#include "TiberModule.h"


using namespace libMesh;



 // ,_incident_angle({0.0})
Tmm::Tmm(const ModelOptions& options) :
  SimulationInterface(options)
{
}

Tmm::Matrix_2by2::Matrix_2by2()
:_m00(0.0),_m01(0.0) ,_m10(0.0) ,_m11(0.0)
{
}

Tmm::Matrix_2by2::Matrix_2by2(double a00, double a01, double a10, double a11)
:_m00(a00),_m01(a01) ,_m10(a10) ,_m11(a11)
{
}

void Tmm::Matrix_2by2::print()
{
  std::cout << "m00 = " << _m00 << "  " << "m01 = "<< _m01 << "\r\n";
  std::cout << "m10 = " << _m10 << "  " << "m11 = "<< _m11 << "\r\n";
  std:cout << "                       " << "\r\n";

}

void Tmm::Matrix_2by2::inv()
{
  complex<double> div (1,0);
  complex<double> m00;
  complex<double> m01;
  complex<double> m10;
  complex<double> m11;

  m00 = _m00;
  m01 = _m01;
  m10 = _m10;
  m11 = _m11;
  div = (m00 *m11 - m01 *m10);
  if (m11 == div)
    _m00 = 1;
  else
  _m00 = m11/div;

  if (m01 == div)
    _m01 = 1;
  else
    _m01 = -m01/div;;

  if (m10 == div)
    _m10 = 1;
  else
  _m10 = -m10/div;

  if (m00 == div)
    _m11 = 1;
  else
  _m11 = m00/div;
}

void Tmm::Matrix_2by2::unit_matrix()
{
  _m00 = 1;
  _m01 = 0;
  _m10 = 0;
  _m11 = 1;
}

Tmm::Matrix_2by2 Tmm::Matrix_2by2::operator*(const Tmm::Matrix_2by2 & old_Matrix_2by2)
{
  Tmm::Matrix_2by2 new_Matrix_2by2;
  new_Matrix_2by2._m00 = (_m00 * old_Matrix_2by2._m00) + (_m01 * old_Matrix_2by2._m10);
  new_Matrix_2by2._m01 = (_m00 * old_Matrix_2by2._m01) + (_m01 * old_Matrix_2by2._m11);
  new_Matrix_2by2._m10 = (_m10 * old_Matrix_2by2._m00) + (_m11 * old_Matrix_2by2._m10);
  new_Matrix_2by2._m11 = (_m10 * old_Matrix_2by2._m01) + (_m11 * old_Matrix_2by2._m11);

  return(new_Matrix_2by2);
}



void Tmm::Matrix_2by2::set(int elm , std::complex<double> a)
{
  switch(elm)
  {
    case 0:
      _m00 = a;
      break;
    case 1:
      _m01 = a;
      break;
    case 2:
      _m10 = a;
      break;
    case 3:
      _m11 = a;
      break;
  }
}

std::complex<double> Tmm::Matrix_2by2::get(int elm)
{
  switch(elm)
  {
    case 0:
      return(_m00);
      break;
    case 1:
      return(_m01);
      break;
    case 2:
      return(_m10);
      break;
    case 3:
      return(_m11);
      break;
    default:
      return(0);
  }
}

std::complex<double> Tmm::cmlx_sqrt(std::complex<double> in)
{
  double real_part=real(in);
  double imag_part=imag(in);

  if (imag_part == 0)
    imag_part = 0.0;
  if (real_part == 0)
    real_part = 0.0;
  std::complex<double> out(real_part,imag_part);
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



Tmm::Matrix_2by2 Tmm::get_M(double n_real,double n_imag,double lenght,double lambda, double kr, double phase)
{
  Tmm::Matrix_2by2 new_Matrix_2by2;
  complex<double> j (0,1);
  complex<double> ki ((2*M_PI*n_real)/lambda , (2*M_PI*n_imag)/lambda);
  complex<double> kz;
  //complex<double> memory;
  kz  = ki * cmlx_sqrt(complex<double> (1.0 - (pow(kr,2)/pow(ki,2))));
  //std::cout<<"kr is :"<<kr<<"   ki :"<<ki<< "   kz :"<<kz<<"    "<<exp(-j * kz)<<std::endl;
  kz *= lenght;
  complex<double> ps (0,phase);

  new_Matrix_2by2.set(0,exp((-j * kz) + ps));
  new_Matrix_2by2.set(1,0);
  new_Matrix_2by2.set(2,0);
  new_Matrix_2by2.set(3,exp((j * kz) - ps));


  return(new_Matrix_2by2);
}

Tmm::Matrix_2by2 Tmm::get_D(double n1_real,double n1_imag,double n2_real,double n2_imag,double kr, double lambda, double mode)
{
  Tmm::Matrix_2by2 new_Matrix_2by2;
  complex<double> j (0,1);
  complex<double> Ni1 (n1_real,n1_imag);
  complex<double> Ni2 (n2_real,n2_imag);
  complex<double> ki1 ((2*M_PI*n1_real)/lambda , (2*M_PI*n1_imag)/lambda);
  complex<double> ki2 ((2*M_PI*n2_real)/lambda , (2*M_PI*n2_imag)/lambda);
  ki1  = ki1 * cmlx_sqrt(complex<double> (1.0 - (pow(kr,2)/pow(ki1,2))));
  ki2  = ki2 * cmlx_sqrt(complex<double> (1.0 - (pow(kr,2)/pow(ki2,2))));
  complex<double> r,t;

  
  if (mode == 0)  //TE Mode
  {
      r   =  (ki1-ki2)/(ki1+ki2);
      t   =  (2.0*ki1)/(ki1+ki2);
  }else           //TM Mode
  {
      r   =   ((-ki1 * Ni2 * Ni2)+(ki2 * Ni1 * Ni1))/((ki1 * Ni2 * Ni2)+(ki2 * Ni1 * Ni1));
      t   =   ( 2 * ki1 * Ni1 * Ni2                )/((ki1 * Ni2 * Ni2)+(ki2 * Ni1 * Ni1));  
  }

  new_Matrix_2by2.set(0,1.0/(t));
  new_Matrix_2by2.set(1,(r)/(t));
  new_Matrix_2by2.set(2,(r)/(t));
  new_Matrix_2by2.set(3,1.0/(t));

  return(new_Matrix_2by2);
}

Tmm::Matrix_2by2 Tmm::Determinal_Matrix (Tmm::Matrix_2by2 MAT)
{
  complex<double> ratio;
  Tmm::Matrix_2by2 MAT_DET;
  ratio = 1/(MAT.get(0)*MAT.get(3) - MAT.get(1)*MAT.get(2));
  MAT_DET.set(0, ratio * MAT.get(3));
  MAT_DET.set(1,-ratio * MAT.get(1));
  MAT_DET.set(2,-ratio * MAT.get(2));
  MAT_DET.set(3, ratio * MAT.get(0));
  return (MAT_DET);

}

void Tmm::dipole_source(double & A_P, double & A_N, double Mode, double Oraintation, double lambda, double cos_phi_inter)
{

  double sin_phi_inter;
  sin_phi_inter = real( cmlx_sqrt(complex<double> (1.0 - pow(cos_phi_inter,2))));


  if (Mode == 0 )    //TE Mode
  {
    if (Oraintation == 0)  // Horizontal 
    {
      A_P = -sqrt(3/(16 * M_PI));
      A_N =  sqrt(3/(16 * M_PI));           
    }else   //Vertical
    {
      A_P = -sqrt(3/(16 * M_PI))*sin_phi_inter;
      A_N =  sqrt(3/(16 * M_PI))*sin_phi_inter; 
    }
  }else     //TM Mode
  {
    if (Oraintation == 0) // Horizontal 
    {
      A_P = -sqrt(3/(16 * M_PI))*cos_phi_inter;
      A_N =  sqrt(3/(16 * M_PI))*cos_phi_inter;             
    }else   //Vertical
    {
      A_P = -sqrt(3/(16 * M_PI))*sin_phi_inter;
      A_N =  sqrt(3/(16 * M_PI))*sin_phi_inter; 
    }

  }
}

void Tmm::solving_internal_source(vector<complex<double>> & E_int_f, vector<complex<double>> & E_int_b ,const vector<double> & n_real,
 const vector<double> & n_imag,const vector<double> & l, double lambda, double dipole_loc, double kr, double Mode, double A_P, double phase)
{
  Tmm::Matrix_2by2 Es(A_P,0,-A_P,0);
  Tmm::Matrix_2by2 DD(1,0,0,1);
  Tmm::Matrix_2by2 MM(0,0,0,0);
  Tmm::Matrix_2by2 TT_load(0,0,0,0);
  Tmm::Matrix_2by2 T_RIGHT(1,0,0,1);
  Tmm::Matrix_2by2 T_LEFT(1,0,0,1);
  Tmm::Matrix_2by2 T_TOTAL(1,0,0,1);
  Tmm::Matrix_2by2 AA(0,0,0,1);

  //*****************************loop for right side of the dipole******************************
  MM = get_M(n_real[n_real.size()-1],n_imag[n_real.size()-1],l[n_real.size()-1],lambda,kr,phase);
  // std::cout<<n_real[n_real.size()-1]<<" , "<<n_imag[n_real.size()-1]<<" , "<< l[n_real.size()-1]<<" , "<< lambda <<" , "<< kr <<std::endl;
  T_RIGHT = MM;
  // MM.print();
  for (double k = n_real.size()-1; k >dipole_loc; --k)
  {
  DD = get_D(n_real[k-1],n_imag[k-1],n_real[k],n_imag[k],kr,lambda,Mode);
  MM = get_M(n_real[k-1],n_imag[k-1],l[k-1],lambda,kr,phase);
  TT_load = DD * T_RIGHT;
  T_RIGHT = MM* TT_load;
  }
  MM = get_M(n_real[dipole_loc-1],n_imag[dipole_loc-1],l[dipole_loc-1]/2,lambda,kr,phase);
  TT_load = MM * T_RIGHT;
  T_RIGHT = TT_load;
  // T_RIGHT.print();

  //*****************************loop for left side of the dipole******************************
  MM = get_M(n_real[dipole_loc-1],n_imag[dipole_loc-1],l[dipole_loc-1]/2,lambda,kr,phase);
  T_LEFT = MM;
  for (double k = dipole_loc-2; k >0; --k)
  {
    DD = get_D(n_real[k],n_imag[k],n_real[k+1],n_imag[k+1],kr,lambda,Mode);
    MM = get_M(n_real[k],n_imag[k],l[k],lambda,kr,phase);
    TT_load = DD * T_LEFT;
    T_LEFT  = MM * TT_load;
  }
  MM = get_M(n_real[0],n_imag[0],l[0],lambda,kr,phase);
  TT_load = MM * T_LEFT;
  T_LEFT = TT_load;
  // T_LEFT.print();
  //*****************************Solving a pair of linear equations to get E.F. in two ends******************************
  T_LEFT.inv();
  Tmm::Matrix_2by2 Sol(1,0,0,1);
  Sol.set(0,T_RIGHT.get(0));
  Sol.set(1,-T_LEFT.get(1));
  Sol.set(2,T_RIGHT.get(2));
  Sol.set(3,-T_LEFT.get(3));
  Sol.inv();
  TT_load = Sol * Es;                 // Electric Field in two ends
  // std::cout<<"Out side E.F is:"<<std::endl;
  //  TT_load.print();
  //*****************************Defining vectors for Electric Field******************************
  vector<complex<double>> E_int_f_r(n_real.size());
  vector<complex<double>> E_int_b_r(n_real.size());
  vector<complex<double>> E_int_f_l(n_real.size());
  vector<complex<double>> E_int_b_l(n_real.size());
  //*****************************Assigning Electric field of two ends of the device******************************
  E_int_f_r[n_real.size()-1] = TT_load.get(0);
  E_int_b_r[n_real.size()-1] = 0;
  E_int_f_l [0] = 0;
  E_int_b_l [0] = TT_load.get(2);
  //*********using last sections results to calculate electric field for the entire device******************
  Tmm::Matrix_2by2 E_N(1,0,0,0);
  Tmm::Matrix_2by2 E_I(0,0,0,0);
  E_N.set(0,E_int_f_r[n_real.size()-1]);
  E_N.set(2,E_int_b_r[n_real.size()-1]);

  T_RIGHT.unit_matrix();
  DD.unit_matrix();
  //*****************************loop for right side TMM******************************
  for (double k = n_real.size(); k >dipole_loc; --k)
  {
    if (k < n_real.size())
      DD = get_D(n_real[k-1],n_imag[k-1],n_real[k],n_imag[k],kr,lambda,Mode);
    MM = get_M(n_real[k-1],n_imag[k-1],l[k-1],lambda,kr,phase);
    TT_load = DD * T_RIGHT;
    T_RIGHT = MM* TT_load;
    E_I = T_RIGHT * E_N;
    E_int_f_r[k-2] = (E_I.get(0));
    E_int_b_r[k-2] = (E_I.get(2));
  }

  E_N.set(0,E_int_b_l[0]);
  E_N.set(2,E_int_f_l[0]);


  T_LEFT.unit_matrix();
  DD.unit_matrix();
  //*****************************loop for left side TMM******************************
  for (double k = 0 ; k <dipole_loc-1; ++k)
  {
    if (k >  0)
      DD = get_D(n_real[k-1],n_imag[k-1],n_real[k],n_imag[k],kr,lambda,Mode);
    MM = get_M(n_real[k],n_imag[k],l[k],lambda,kr,phase);
    DD.inv();
    TT_load = DD * T_LEFT;
    T_LEFT = MM* TT_load;
    E_I = T_LEFT * E_N;
    E_int_f_l[k+1] = (E_I.get(2));
    E_int_b_l[k+1] = (E_I.get(0));
  }

  //*****************************loops to to sum up electric fields******************************

  for (double nm = 0; nm<dipole_loc-1;nm++)
  {
    E_int_f[nm] = E_int_f_l[nm+1];
    E_int_b[nm] = E_int_b_l[nm+1];
  }
  for (double nm = E_int_f_r.size()-1; nm>dipole_loc-1;nm--)
  {
    E_int_f[nm] = E_int_f_r[nm-1];
    E_int_b[nm] = E_int_b_r[nm-1];
  }
  E_int_f[dipole_loc-1] = (E_int_f[dipole_loc-2]+E_int_f[dipole_loc])/2;
  E_int_b[dipole_loc-1] = (E_int_b[dipole_loc-2]+E_int_b[dipole_loc])/2;
}



vector<double> Tmm::theta_cal(vector<double> n_real , double incident_angle)
{
  vector<double> theta(n_real.size());
  theta[0]=incident_angle;
  for (int k=1; k<n_real.size();k++)
  {
    theta[k]=sin(theta[k-1]*M_PI/180);

    theta[k]=(n_real[k-1]/n_real[k])*theta[k];
 
    theta[k]=asin(theta[k])*180/M_PI;
  
  }
  return theta;
}


vector<double> Tmm::linear_interpolation1 (vector<double> xData, vector<double> yData, vector<double> x_interp)
{
  int size = xData.size();
  vector<double> y_interp;

  for (double x : x_interp)
  {
    int i =0;
    if (x >= xData[size-2])
    {
      i = size -2;
    }
    else
    {
      while (x > xData[i+1]) i++;
    }
    double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];
    double delta = (yR - yL) / (xR - xL);
    double y = yL + delta * (x-xL);
    y_interp.push_back(y);
  }
  return y_interp;
}



Tmm::~Tmm(void)
{
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

  _incident_angle = get_option("incident_angle",0);

  _reflectivity = get_option("back_reflectivity",0);

  _dipole_loc = get_option("dipole_loc",0);


  _coh_mod = get_option("coh_mode",0);



  _polarization = get_option("polarization","");
  _orientation  = get_option("orientation" ,"");
  
  if (_polarization.empty())
    _polarization_vec.push_back(0);

  if (cmp_string(_polarization,"TEM"))
  {
    _polarization_vec.push_back(0);
    _polarization_vec.push_back(1);
  }
  else
  {
  if (cmp_string(_polarization,"TE"))
      _polarization_vec.push_back(0);
    else
      if (cmp_string(_polarization,"TM"))
        _polarization_vec.push_back(1);
      else
        Messages::warning("Unknown polarization.");
  }



  if (_orientation.empty())
    _oraintation_vec.push_back(0);

  if (cmp_string(_orientation,"VH"))
  {
    _oraintation_vec.push_back(0);
    _oraintation_vec.push_back(1);
  }
  else
  {
    if (cmp_string(_orientation,"H"))
      _oraintation_vec.push_back(0);
    else
      if (cmp_string(_orientation,"V"))
        _oraintation_vec.push_back(1);
      else
        Messages::warning("Unknown orientation.");
  }


  _steps = get_option("wave_number_steps",0);
  get_option("wave_number_ratio",_ratio);





  get_option("dipole_power",_dipole_power);
  get_option("dipole_coordinate",_dipole_coordinate );

  if (_dipole_coordinate.size() != _dipole_power.size() )
     if (_dipole_power.size() == 1 )
     {
      _dipole_power.resize(_dipole_coordinate.size());
      for (int nm = 1 ; nm < _dipole_coordinate.size(); nm++ )
        _dipole_power[nm] = _dipole_power[0];
     }
     else
      Messages::warning("The size of dipole_power and dipole_coordinate vectors don't match!");

  get_option("wavelengths",_wavelength_vector );
  if (_wavelength_vector.empty())
  {
      _up_lambda = get_option("wavelength_uper_lim", 0 );
      if (_up_lambda == 0)
      {
        Messages::warning("You did not provide any up_lambda for TMM.");
      }
      _down_lambda = get_option("wavelength_lower_lim", 0);
      if (_down_lambda == 0)
      {
        Messages::warning("You did not provide any down_lambda for TMM.");
      }
      _wavelength_steps = get_option("wavelength_steps",1);
  }


  Database db;
  ifstream is;
  db.set_material("Sun1p5am", get_option("illumination_spectrum", ""));
  is.open(db.get_data_file().c_str());
  if (is.fail() || !is.good())
    throw InitFailedException("Cannot read spectrum "
        "from file " + db.get_data_file());

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


void
Tmm::do_setup_solution_variables(void)
{
  // we declare our solution variables
  declare_solution(GenerationRate, REAL, CELL, "cm^3");
  declare_solution(Intensity, REAL, CELL, "mW/cm^3nm^-1");
  declare_solution(External_Source_ElectricField, REAL, CELL, "V/cm");

  declare_solution(Internal_Source_ElectricField, REAL, CELL, "V/cm");
  declare_solution(Internal_Poynting, REAL, CELL, "a.u.");
  declare_solution(Internal_Power, REAL, CELL, "a.u.");
  declare_solution(Internal_Absorption, REAL, CELL, "a.u.");
  declare_solution(Internal_Intensity, REAL, CELL, "a.u.");

  declare_solution(Transmission, REAL, GLOBAL, "1");
  declare_solution(Reflection, REAL, GLOBAL, "1");
  declare_solution(Absorption, REAL, GLOBAL, "1");

}


void
Tmm::do_solve(void)
{
  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  NumericVector<libMesh::Number>& solution = system.get_local_solution_vector();
  solution.close();
  solution.zero();


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

  DofMap& dof_map =  system.get_dof_map();
  double mesh_units = get_mesh_units();
  vector<unsigned int> dof_indices;
  vector<double> intensity_integral;
  vector<double> Electric_Field_integral;
  vector<double> generation_rate_integral;
  vector<double> gen;
  vector<double> lambda_interp;
  vector<double> sun_interp;
  vector<double> Regions;

  double external_source_simulation=0;
  double internal_source_simulation=0;
  vector<complex<double>> Ki;

  if (_wavelength_vector.empty())
    for (double i = _down_lambda; i<= _up_lambda; i += _wavelength_steps)
      lambda_interp.push_back(i);
  else
    for (double i = 0; i< _wavelength_vector.size(); ++i)
      lambda_interp.push_back(_wavelength_vector[i]);

  sun_interp = Tmm::linear_interpolation1(_lambda,_spectrum,lambda_interp);


  for (unsigned int i = 0; i < lambda_interp.size(); ++i)   //loop over wavelength
  {
    ostringstream os;
    os << "----------------------------------------"<<"\n" << "solving condition :  " ;
    os << "Lambda is : " << lambda_interp[i] << "nm" << "\n" ;
    Messages::info(os.str());

    double lambda = lambda_interp[i];
    _Wavelength.push_back(lambda);



    const unsigned int uvar = system.variable_number("G");


    // TODO this will not work if the 1D mesh is distributed. In that case, MPI calls could be used
    // to gather pieces from all processes
    MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

    double E = 1.0;
    double incoming_angle=_incident_angle;

    double c0 = 2.998e8 * 1e9;
    double e0 = 8.85e-12;
    double radiation = sun_interp[i] / 1e3;
    double Esun2 = 2 * radiation / (c0 * 1e-9 * e0);
    double w = 2 * M_PI * c0 / lambda;
    double plank_const = 1.055e-34;




    // TODO reserve space
    vector<double> n_real;
    vector<double> n_imag;
    vector<double> l_length;
    vector<double> l;
    vector<double> Incoh;
    vector<double> elem_coordinate;;


    vector<double> n_real_init;
    vector<double> n_imag_init;
    vector<double> l_length_init;
    vector<double> l_init;
    vector<double> Incoh_init;

    vector<std::string> Names;
    


    vector<int> BC_check;
    vector<int> BC_check_init;
    Tmm::Matrix_2by2 Boundry_condition(1,0,0,1);

    std::string direction;


    //**********************************************************************************************

    for ( ; el != end_el ; ++el)
    {


      const Elem* elem = *el;
      dof_map.dof_indices(elem, dof_indices, uvar);
      const unsigned int n_dofs = dof_indices.size();
      TmmBulkModel& mod = *get_bulk_model<TmmBulkModel>(elem);
      mod.reinit(elem);
      l_length.push_back(dof_indices[0]);
      libMesh::Point pp =elem->centroid();
      elem_coordinate.push_back(pp(0));

      if ( i == 0)
      {
       ID sub_id = elem->subdomain_id();
       std::string name = get_environment().get_device().get_region_name(sub_id);
       Names.push_back(name);
      }

      libMesh::Complex nk = mod.get_refractive_index(lambda);
      Incoh_init.push_back(mod.get_coherent_index());
      n_real_init.push_back(real(nk));
      n_imag_init.push_back(imag(nk));
      l_init.push_back(elem->volume()*(mesh_units/1e-9));
       for (unsigned int s = 0; s < elem->n_sides(); s++)
       {
         TmmBoundaryModel* mod_int =
           get_interface_model<TmmBoundaryModel>(elem, s);
         if (mod_int != NULL)
         {
           /*
           if (internal_source_simulation==0 && mod_int->read_type() == "Dipole Source"){
             internal_source_simulation = 1;
             kr_ratio = mod_int->get_kr();
             steps = mod_int->get_steps();
             if (steps == 0)
               steps = 1;
             for (double nm =0; nm<l_init.size();nm++)
               _dipole_loc += l_init[nm];
           }
           */
           if (mod_int->read_type() == "Mirror"){
             mod_int->Calculate_M_Matrix();
             BC_check_init.push_back(1);
             Boundry_condition.set(0,mod_int->get_element(0));
             Boundry_condition.set(1,mod_int->get_element(1));
             Boundry_condition.set(2,mod_int->get_element(2));
             Boundry_condition.set(3,mod_int->get_element(3));
           }
           if (mod_int->read_type() == "Incident Wave"){
             external_source_simulation = 1;
             BC_check_init.push_back(0);
             if ( s % 2 == 0)
             {
               direction = "top to down propagation";
             }
             else
             {
               direction = "down to top propagation";
             }
           }
         }
           if (mod_int == NULL)
           {
             //std::cout<<"no boundry condition found" <<std::endl;
             BC_check_init.push_back(0);
           }
       }
    }

    if (direction == "down to top propagation")
    {
      for (int uu=BC_check_init.size()-2; uu>=0;uu -= 2)
        BC_check.push_back(BC_check_init[uu]);
      for (int uu=n_real_init.size()-1; uu>=0;--uu)
      {
        n_real.push_back(n_real_init[uu]);
        n_imag.push_back(n_imag_init[uu]);
        l.push_back(l_init[uu]);
        Incoh.push_back(Incoh_init[uu]);
      }

    }else
    {
      for (int uu=1; uu<BC_check_init.size();uu += 2)
        BC_check.push_back(BC_check_init[uu]);
      for (int uu=0; uu<n_real_init.size();++uu)
      {
        n_real.push_back(n_real_init[uu]);
        n_imag.push_back(n_imag_init[uu]);
        l.push_back(l_init[uu]);
        Incoh.push_back(Incoh_init[uu]);
      }
    }

    if ( i == 0)
    {
      for (double lm = 0; lm < Names.size();lm++)
        if (lm ==0 || Names[lm] != Names[lm-1] )
        {
          _regions_name.push_back(Names[lm]);
          Regions.push_back(lm);
        }
      Regions.push_back(n_real.size());
    }
   

    //****************************************************************************
    //*****************************External_source_simulation******************************
    if (external_source_simulation)
    {
      
      vector<double> theta(n_real.size());
      theta=Tmm::theta_cal(n_real,incoming_angle);
      double rnd = 1;
      double phase_step;
      vector<double> Generation_rate_avg(n_real.size());
      vector<double> Intensity_avg(n_real.size());
      vector<double> Electric_Field_avg(n_real.size());
      for(int nm=0; nm<Intensity_avg.size();++nm)
      {
        Intensity_avg[nm] = 0;
        Generation_rate_avg [nm] = 0;
        Electric_Field_avg [nm] = 0;
      }

      double avg_reflection = 0;
      double avg_transmission = 0;
      double avg_Absorption = 0;


      double coh_mod = _coh_mod;

      if(coh_mod)
      {
        for (int nm = 0; nm<Incoh.size(); ++nm)
          if(Incoh[nm]==1)
            rnd = _coh_mod;
      }
      else
      {
        for (int nm = 0; nm<Incoh.size(); ++nm)
          if(Incoh[nm]==1)
            rnd = 5;
      }

      phase_step = 2 * M_PI / rnd;
      double N0 = 0;
      srand(time(NULL));
      for (int iter = 0; iter <rnd; ++iter )    // loop over added phases
      {
        
        double phase;
        double r =_reflectivity;
        //***********************************************************************
        //**********loop added for highly absorbing devices*******************
        Tmm::Matrix_2by2 Unit(1,0,0,1);
        Tmm::Matrix_2by2 DD(1,0,0,1);
        Tmm::Matrix_2by2 MM(0,0,0,0);
        Tmm::Matrix_2by2 TT_load(0,0,0,0);
        Tmm::Matrix_2by2 TT(1,0,0,1);
        for (int k = 0 ;k <= n_real.size()-2 ; ++k)
        {
          DD = get_D(n_real[k],n_imag[k],n_real[k+1],n_imag[k+1],0,lambda,0);

          if(Incoh[k]==1 && Incoh[k+1]==0)
          {
            MM = get_M(n_real[k],n_imag[k],l[k],lambda,0,phase);
          }
          else
          {
            MM = get_M(n_real[k],n_imag[k],l[k],lambda,0,0);
           }
           TT_load = MM * DD;
           if (real(abs(TT.get(0))) > 1e150)
           {
             TT_load = Unit;
             TT = TT * TT_load ;
             if (N0 == 0)
               N0 = n_real.size() - k;
            }else
              TT = TT * TT_load ;
          }
        //***********************************************************************
        //**********************Main TMM loop************************************
          Tmm::Matrix_2by2 T(1,0,0,1);
          Tmm::Matrix_2by2 E_N(1,0,0,0);
          Tmm::Matrix_2by2 E_I(0,0,0,0);
          vector<complex<double>> E_F;
          vector<complex<double>> E_B;

          if (N0 != 0)
               for (int nm = 0 ; nm< N0; nm++)
               {
                 E_F.push_back(0);
                 E_B.push_back(0);
               }
           E_F.push_back(1);
           E_B.push_back(0);
           //T = get_M(n_real[n_real.size()-N0-2],n_imag[n_real.size()-N0-2],l[n_real.size()-N0-2],lambda,theta[n_real.size()-N0-2],0);

           for (int k = n_real.size()-N0-2 ;k >= 0 ; --k)
           {
             DD = get_D(n_real[k],-n_imag[k],n_real[k+1],-n_imag[k+1],0,lambda,0);
             if(Incoh[k]==1 && Incoh[k+1]==0)
             {
              if (coh_mod)
              {
                
                phase = 2 * M_PI *(double) rand()/RAND_MAX;
              }
              else
              {
                phase = phase_step * iter;
              }
              MM = get_M(n_real[k],n_imag[k],l[k],lambda,0,phase);
             }
             else
             {
              MM = get_M(n_real[k],n_imag[k],l[k],lambda,0,0);
             }
             TT_load = DD * T;
             T = MM* TT_load;
             E_I = T * E_N  ;
             E_F.push_back(E_I.get(0)) ;
             E_B.push_back(E_I.get(2));
           }
           //****************************************************************************
           //***************normalizing electric field matrix****************************
           vector<complex<double>> E_F_NORM(n_real.size());
           vector<complex<double>> E_B_NORM(n_real.size());

           for(double nm=0; nm < n_real.size() ; ++nm)
           {
             E_F_NORM[nm] = E_F[E_F.size()-1-nm]/E_F[E_F.size()-1];
             E_B_NORM[nm] = E_B[E_F.size()-1-nm]/E_F[E_F.size()-1];
           }
        //***********************************************************************
        //**********reflection and transmission calculation**********************

        complex<double> Reflection,Transmission;
        Reflection = pow(abs(T.get(2)/T.get(0)),2);
        avg_reflection = avg_reflection + real(Reflection) / rnd;

        complex<double> nc_first (n_real[0],n_imag[0]);
        complex<double> nc_last (n_real[n_real.size()-1],n_imag[n_imag.size()-1]);
        complex<double> ratio_complex;
        ratio_complex = ((nc_last)*cos(theta[theta.size()-1]*M_PI/180))/(nc_first*cos(theta[0]*M_PI/180));
        Transmission = ratio_complex*pow(abs(1.0/T.get(0)),2);
        // std::cout<<"T is : "<<Transmission<<std::endl;
        avg_transmission = avg_transmission + real(Transmission) / rnd;
        avg_Absorption = avg_Absorption + (1-real(Reflection)-real(Transmission)) / rnd;

        //****************************************************************************
        //***************Calculating average values over random phases****************

        vector<complex<double>> Etot(n_real.size());
        vector<complex<double>> Intensity;
        vector<complex<double>> Generation_rate;
        vector<double> Generation_rate_real;

        for (int nm=0 ; nm < n_real.size() ; ++nm)
        {
          Etot[nm] = E_F_NORM[nm] + E_B_NORM[nm];
          Electric_Field_avg[nm] =  Electric_Field_avg[nm] + real(Etot[nm])/rnd;
          Intensity.push_back(0.5 * c0 * 1e-9 * e0 * n_real[nm] * Esun2* pow(abs(Etot[nm]), 2));
          Intensity_avg[nm] = Intensity_avg[nm] + real(Intensity[nm])/rnd;
          Generation_rate.push_back(1 / (plank_const* w) * (4 * M_PI * n_imag[nm] * 1e7/(lambda)) * real(Intensity[nm])/ 1e4 );
          Generation_rate_real.push_back(real(Generation_rate[nm]));
          Generation_rate_avg[nm] = Generation_rate_avg[nm] + Generation_rate_real[nm]/rnd;
        }
      }
      // std::cout<<"Totla T  is : "<<avg_transmission<<std::endl;
      vector<double> sumation;
      for(double re =1; re < Regions.size(); re++)
      {
        double load=0;
        for (double el=Regions[re-1]; el < Regions[re]; el++)
        {
          load += Generation_rate_avg[el];
        }
        sumation.push_back(load);
      }
      double load=0;
      for (double re = 0; re<sumation.size();re++)
        load += sumation[re];
      for (double re = 0; re<sumation.size();re++)
        sumation[re] = (sumation[re]/load) * avg_Absorption;

      _Generation_regions.push_back(sumation);

      _Electric_Field_External.push_back(Electric_Field_avg);
      
      double scale = 1e-6;
      _Reflection.push_back((int)(real(avg_reflection) / scale) * scale);
      _Transmission.push_back((int)(real(avg_transmission) / scale) * scale);
      _Absorption.push_back((int)(avg_Absorption / scale) * scale);

      //****************************************************************************
      //*******************Calculating integral over wavelengths********************
      if (i == 0) // first wavelength
      {

        Electric_Field_integral.resize(n_real.size());
        Electric_Field_integral = Electric_Field_avg;

        intensity_integral.resize(n_real.size());
        intensity_integral = Intensity_avg;

        generation_rate_integral.resize(n_real.size());
        generation_rate_integral = Generation_rate_avg;

      }else // middle wavelengths
      {
        for (int nm =0; nm< n_real.size();nm++)
        {
          generation_rate_integral[nm] +=  Generation_rate_avg[nm];
          intensity_integral[nm] += Intensity_avg[nm];
          Electric_Field_integral[nm] += Electric_Field_avg[nm];
        }
      }
      if( i == lambda_interp.size()-1) //last wavelength
      {
         i++ ;
        for (int nm =0; nm< n_real.size();nm++)
        {
          if (i !=1)
          generation_rate_integral[nm] *= (lambda_interp[lambda_interp.size()-1]-lambda_interp[0])/(i);

          intensity_integral[nm] *= 1.0/i;
          Electric_Field_integral[nm] *= 1.0/i;
        }
        for(double nm=0; nm < l_length.size() ; ++nm){
          _Generation_rate.push_back(generation_rate_integral[nm]);
          _Intensity.push_back(intensity_integral[nm]);
          _External_Source_ElectricField.push_back(Electric_Field_integral[nm]);
        }
      }
    }
    //*****************************End of External_source_simulation******************************

    //****************************************************************************
    //*****************************Internal_source_simulation******************************


    if (!_dipole_coordinate.empty()) // if any dipole is defined, solve it.
      internal_source_simulation = 1;

    if (internal_source_simulation)
    {
      for (double dipole_num=0 ;dipole_num<_dipole_coordinate.size();++ dipole_num)
        for (double elem = 1 ; elem < elem_coordinate.size(); elem++)
          if (_dipole_coordinate[dipole_num] >= elem_coordinate[elem-1] && _dipole_coordinate[dipole_num] < elem_coordinate[elem])
            _dipole_coordinate[dipole_num] = elem-1;


      vector<double> kr_vec;
 
      double ks = 2*M_PI*n_real[_dipole_coordinate[0]]/lambda;
      if (_steps != 0)
      {
        kr_vec.resize(_steps);
        kr_vec[0] = 0;
        for (double u = 1; u < _steps; ++u)
        {
          kr_vec[u] = kr_vec[u-1] + ks*_ratio[0]/_steps;
        }
      }else
      {
        _steps = 1;
        kr_vec.resize(_steps);
        kr_vec[0] = ks*_ratio[0];
      }

      double A_P ;
      double A_N ; 
      
      double rnd = 1;
      double phase_step;
      
      for (int nm = 0; nm<Incoh.size(); ++nm)
        if(Incoh[nm]==1)
          rnd = 5;
      
      phase_step = 2 * M_PI / rnd;
      _Internal_Source_ElectricField.resize(n_real.size());
      _Internal_Poynting.resize(n_real.size());
      _Internal_Power.resize(n_real.size());
      _Internal_Absorption.resize(n_real.size());
      _Internal_Intensity.resize(n_real.size());
      _Fraction_ratio.resize(_steps);
      _kr.resize(_steps);
      _angle.resize(_steps);
      _Poynting_front.resize(_steps);
      _Poynting_back.resize(_steps);



      for (double polarization_num = 0 ; polarization_num < _polarization_vec.size() ; polarization_num ++)
      {
        for (double oraintation_num = 0 ; oraintation_num < _oraintation_vec.size() ; oraintation_num ++)
        { 
          if (_polarization_vec[polarization_num] == 0 )    //TE Mode
                if (_oraintation_vec[oraintation_num] == 0)  // Horizontal 
                  std::cout<<"TE Mode and Horizontal"<<std::endl;         
                else   //Vertical
                  std::cout<<"TE Mode and Vertical"<<std::endl;  
          else     //TM Mode
                if (_oraintation_vec[oraintation_num] == 0) // Horizontal 
                  std::cout<<"TM Mode and Horizontal"<<std::endl;            
                else   //Vertical
                  std::cout<<"TM Mode and Vertical"<<std::endl; 


        //*****************************loop for number of dipoles******************************
          for (double dipole_num=0 ;dipole_num<_dipole_coordinate.size();++ dipole_num)
          {
            double dipole_loc = _dipole_coordinate[dipole_num];
            std::cout<<"Solving for internal source --> dipole coordinate is "<<elem_coordinate[dipole_loc]<<" internal power is : "<<_dipole_power[dipole_num] <<std::endl;
            complex<double> Ns (n_real[dipole_loc], n_imag[dipole_loc]);
           //*****************************loop for radial wave numbers******************************
            double cnt =0;
            double kr;
            for (double kr_num = 0; kr_num < _steps; ++kr_num )
            {
              double Mode = _polarization_vec[polarization_num];    //TE Mode
              double Oraintation = _oraintation_vec[oraintation_num]; //Vertical Mode
              kr = kr_vec[kr_num];
              // std::cout<<"kr is : "<<kr << "  , cnt is "<< cnt <<std::endl;

              if (dipole_num == 0)
              {
                _kr[cnt] = kr;
                _Fraction_ratio[cnt] = abs(kr/ks);
              }
                
              complex<double> ki_s    ((2*M_PI*n_real[dipole_loc]) / lambda , (2*M_PI*n_imag[dipole_loc])  / lambda); 
              double cos_phi_inter;
              cos_phi_inter  = real(cmlx_sqrt(complex<double> (1.0 - (pow(kr,2)/pow(ki_s ,2))))) ;
              if (cos_phi_inter > 1)
                cos_phi_inter = 1;
              _angle [cnt] = acos(cos_phi_inter );
              
              dipole_source(A_P, A_N, Mode, Oraintation, lambda, cos_phi_inter);

              

            //*****************************calculating electric field of two ends of the device******************************
              vector<complex<double>> E_int_f(n_real.size());
              vector<complex<double>> E_int_b(n_real.size());
              vector<complex<double>> E_int(n_real.size());
              for (int iter = 0; iter <rnd; ++iter )    // loop over added phases
              {
                solving_internal_source(E_int_f, E_int_b ,n_real,n_imag, l, lambda, dipole_loc,  kr, Mode, A_P,phase_step * iter);
                for (double nm =0; nm < E_int_f.size(); ++nm)
                  E_int[nm] =  E_int[nm] + (E_int_f[nm]+E_int_b[nm])/rnd;
              }
                

              double coefficient;
              coefficient = _dipole_power[dipole_num]/_steps/lambda_interp.size()/_oraintation_vec.size()/_polarization_vec.size();
              for (double nm =0; nm < E_int_f.size(); ++nm)
              {
                // E_int_f[nm] /= rnd;
                // E_int_b[nm] /= rnd;
                // E_int[nm] = E_int_f[nm]+E_int_b[nm];
                _Internal_Source_ElectricField[nm] += real(E_int[nm])         * coefficient;
                _Internal_Intensity[nm] += n_real[nm] * abs(pow(E_int[nm],2)) * coefficient;
              }


              if ( kr_num == 0 && dipole_num == 0 && oraintation_num==0 && polarization_num ==0)
              {
                std::vector<double> buffer(E_int.size());
                for (double cnter =0 ; cnter < E_int.size();cnter++)
                  buffer[cnter] = real(E_int[cnter]);
                _Electric_Field_Internal.push_back(buffer);
              }
       

                //**************a loop to calculate Poynting Vector from electric field******************************
                vector<double> poynting(n_real.size());
                for (double nm = 0; nm<E_int.size();nm++)
                  {
                      complex<double> ki ((2*M_PI*n_real[nm])/lambda , (2*M_PI*n_imag[nm])/lambda);
                      complex<double> kzi (1,0);
                      complex<double> j (0,1);

                      kzi  = ki * cmlx_sqrt(complex<double> (1.0 - (pow(kr,2)/pow(ki,2))));
                      complex<double> cos_phi;
                      cos_phi= kzi/ki;
                      complex<double> Ni (n_real[nm], n_imag[nm]);
                      if (Mode == 0 )  //TE Mode
                        poynting[nm] = real(Ni*cos_phi      *conj(E_int_f[nm]+E_int_b[nm])*(E_int_f[nm]-E_int_b[nm]));
                      else    //TM Mode
                        poynting[nm] = real(Ni*conj(cos_phi)*conj(E_int_f[nm]+E_int_b[nm])*(E_int_f[nm]-E_int_b[nm]));
                  }

                //*****************************Assigning Output Power of two ends of the device ******************************

                double initial_poynting;
                _Poynting_back [cnt] += abs(poynting[0]                ) * coefficient ;
                _Poynting_front[cnt] +=  abs(poynting[n_real.size()-1] ) * coefficient ;

                //*************a loop to take integral of the absorption and Poynting vector over wavelengths, polarization and orientation******

                for (double nm = 0; nm < poynting.size(); ++nm)
                {
                  if (nm < dipole_loc-1 || nm > dipole_loc+1 )
                    _Internal_Absorption[nm] += real(1 / (plank_const* w) * (4 * M_PI * n_imag[nm] * 1e7/(lambda)) * (abs(poynting[nm])/ 1e4));
                  _Internal_Poynting  [nm] += poynting[nm] * coefficient;
                }
                

                //*****************a loop to calculate Power by taking integral over kr ******************************
                complex<double> Ns ( n_real[dipole_loc] , n_imag[dipole_loc] );
                complex<double> ki ((2*M_PI*n_real[dipole_loc])/lambda , (2*M_PI*n_imag[dipole_loc])/lambda);
                complex<double> kzs;
                kzs  = ki * cmlx_sqrt(complex<double> (1.0 - (pow(kr,2)/pow(ki,2))));
                double dkr = ks*_ratio[0]/_steps;
                for (double nm = 0; nm < n_real.size(); ++nm)
                  _Internal_Power[nm] += real( (2*M_PI) * (poynting[nm]* coefficient) * kr * dkr / Ns / pow(kzs,2) ) ;
                ++cnt;

            } // end of radial wave numbers loop

          } // end of number of dipoles loop

         } // end of number of orain loop

     }// end of number of modes loop


    }// end of dipole simulation


  }// end of loop of wavelength
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
  
  string outdir = get_output_directory();


  if (!_Generation_regions.empty())
  {


    string filename(outdir + "/" + get_output_filename() + "_Absorption_Regions.dat");

    ofstream file;
    file.open(filename.c_str());

    if (file.good())
    {
      file << "# " << get_type() << " Absorption_Regions (" << get_name() << ")\n";
      file << "# " << 0 << " WaveLength [nm]" << "\n";
      for (double i = 0; i < _regions_name.size(); ++i)
      {
          file <<"# " << i+1 << " Region => "<< _regions_name[i] << "\n"; 
      }

      file << "# " << "WaveLength ";
      for (double i = 0; i < _regions_name.size(); ++i)
      {
          file << _regions_name[i] << " "; 
      }
      file << "\n";

      for (unsigned int i = 0; i < _Wavelength.size(); i++)
      {
          file << _Wavelength[i]   << " " ;
          for (unsigned int j =0; j < _regions_name.size() ; j++)
            file << _Generation_regions[i][j]   << " " ;
          file << "\n";
      }
    }

      
    file.close();
  }

  if (!_Electric_Field_External.empty())
  {
    string filename(outdir + "/" + get_output_filename() + "_electric_field_external.dat");

    ofstream file;
    file.open(filename.c_str());

    if (file.good())
    {
      long int step=1;
      std::vector<double> wl_vect;
      if (_Wavelength.size()<10)
      {
        wl_vect.resize(_Wavelength.size());
        for (double i =0; i< _Wavelength.size() ; i++)
          wl_vect[i] = _Wavelength[i];
      }
      else
      {
        wl_vect.resize(10);
        step = (_Wavelength.size()) / 9;
        for (double i = 0; i< 10; i++)
          wl_vect[i] = _Wavelength[step * i];

      }

      // header
      file << "# " << get_type() << " TMMM (" << get_name() << ")\n";
      for (double i = 0; i < wl_vect.size(); ++i)
      {
          file <<"# " << i << " WaveLength = "<< wl_vect[i] << "[nm]" << "\n"; 
      }
      file << "# " << "x ";
      for (double i = 0; i < wl_vect.size(); ++i)
      {
          file << wl_vect[i] << "[nm] "; 
      }
      file << "\n"; 

      for (double j = 0; j<_External_Source_ElectricField.size(); ++j)
      {
        file << j   << " ";
        for (double i = 0; i < wl_vect.size(); ++i)
        {
            file << _Electric_Field_External[i*step][j] << " ";
        }
        file << "\n";
      }
    }
    file.close();
  }

if (!_Electric_Field_Internal.empty())
  {
    string filename(outdir + "/" + get_output_filename() + "_electric_field_inernal.dat");

    ofstream file;
    file.open(filename.c_str());

    if (file.good())
    {
      long int step=1;
      std::vector<double> wl_vect;
      if (_Wavelength.size()<10)
      {
        wl_vect.resize(_Wavelength.size());
        for (double i =0; i< _Wavelength.size() ; i++)
          wl_vect[i] = _Wavelength[i];
      }
      else
      {
        wl_vect.resize(10);
        step = (_Wavelength.size()) / 9;
        for (double i = 0; i< 10; i++)
          wl_vect[i] = _Wavelength[step * i];

      }
      // header
      file << "# " << get_type() << " TMMM (" << get_name() << ")\n";
      for (double i = 0; i < wl_vect.size(); ++i)
      {
          file <<"# " << i << " WaveLength = "<< wl_vect[i] << "[nm]" << "\n"; 
      }
      file << "# " << "x ";
      for (double i = 0; i < wl_vect.size(); ++i)
      {
          file << wl_vect[i] << "[nm] "; 
      }
      file << "\n"; 

      for (double j = 0; j<_Internal_Source_ElectricField.size(); ++j)
      {
        file << j   << " ";
        for (double i = 0; i < wl_vect.size(); ++i)
        {
            file << _Electric_Field_Internal[i*step][j] << " ";
        }
        file << "\n";
      }
    }
    file.close();
  }


  if (!_Electric_Field_External.empty())
  {
    string filename(outdir + "/" + get_output_filename() + ".dat");

    ofstream file;
    file.open(filename.c_str());

    if (file.good())
    {
      // header
      file << "# " << get_type() << " TMMM (" << get_name() << ")\n";

      file << "# WaveLength[nm] " <<  " Transmission[1] "
          << " Reflection[1] "
           << " Absorption[1] " <<  "\n";

      for (unsigned int i = 0; i < _Transmission.size(); i++)
      {
          file << _Wavelength[i]   << " " 
               << _Transmission[i] << " "
               << _Reflection[i]   << " "
               << _Absorption[i]   << "\n";
      }
    }
    file.close();
  }


  if (!_Fraction_ratio.empty())
  {
    string polar_file(outdir + "/" + get_output_filename() + "_polar.dat");
    ofstream polar;
    polar.open(polar_file.c_str());
    if (polar.good())
      polar << '#' << get_type() << " Radiation results with respect to radial wave number (" << get_name() << ")\n"
            << "# 0 ratio[1] "             << "\n"
            << "# 1 kr value[1]"           << "\n"
            << "# 2 Emission angle [deg]"  << "\n"
            << "# 3 Power flow back"       << "\n"
            << "#  Powerflow front"        << "\n"
            << "# ratio  kr  Emission_angle  Power_back  Power_front"<< "\n";


    for (double i = 0; i < _Poynting_back.size(); ++i)
      polar << _Fraction_ratio[i] << "  "
            << _kr[i]             << "  "
            << _angle[i]          << "  "
            << _Poynting_back[i]  << "  "
            << _Poynting_front[i] << "\n";

    polar.close();
  }
}

void
Tmm::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& solutions,
    const std::vector<Point>& p)
{
  unsigned int np = p.size();

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();
  const NumericVector<libMesh::Number>& solution = system.get_solution_vector();
  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system.get_dof_map();

  const unsigned int u_var = system.variable_number("G");

  FEType fe_type = system.variable_type(u_var);
  UniquePtr<FEBase> fe(build_finite_element(dim, fe_type));

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

  if (!_Internal_Intensity.empty())
  {
    solutions[Internal_Intensity][0]= _Internal_Intensity[dof_indices[0]];
  }

  if (!_Internal_Absorption.empty())
  {
    solutions[Internal_Absorption][0]= _Internal_Absorption[dof_indices[0]];
  }

  if (!_Internal_Poynting.empty())
  {
    solutions[Internal_Poynting][0]= _Internal_Poynting[dof_indices[0]];
  }

  if (!_Internal_Power.empty())
  {
    solutions[Internal_Power][0]= _Internal_Power[dof_indices[0]];
  }

  if (!_Internal_Source_ElectricField.empty())
  {
    solutions[Internal_Source_ElectricField][0]= _Internal_Source_ElectricField[dof_indices[0]];
  }

  if (!_External_Source_ElectricField.empty())
  {
    solutions[External_Source_ElectricField][0]= _External_Source_ElectricField[dof_indices[0]];
  }

  if (!_Intensity.empty())
  {
    solutions[Intensity][0]= _Intensity[dof_indices[0]];
  }

  if (!_Generation_rate.empty())
  {
    solutions[GenerationRate][0] =_Generation_rate[dof_indices[0]];
  }

}

