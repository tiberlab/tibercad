// $Id$

#include "Tmm.h"
#include "TiberLinearSystem.h"
#include "TmmBulkModel.h"
#include "TmmBoundaryModel.h"
#include "Messages.h"
#include "Database.h"


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
  complex<double> load(1,0);
  div = (_m00 *_m11 - _m01 *_m10);
  load = _m00;
  _m00 = _m11/div;
  _m01 = -_m01/div;
  _m10 = -_m10/div;
  _m11 = load/div;
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


Tmm::Matrix_2by2 Tmm::get_M(double n_real,double n_imag,double lenght,double lambda, double kr, double phase)
{
  Tmm::Matrix_2by2 new_Matrix_2by2;
  complex<double> j (0,1);
  complex<double> ki ((2*M_PI*n_real*lenght)/lambda , (2*M_PI*n_imag*lenght)/lambda);
  //std::cout<<" kr is " << kr <<"  ki is " << ki;
  ki = ki * sqrt(1-(pow(kr,2)/pow(ki,2)));
  //bi=bi* cos(kr*M_PI/180);
 //std::cout<<"  kzi is " <<ki<<" and " <<exp(-j * ki)<<std::endl;
  complex<double> ps (0,phase);


  new_Matrix_2by2.set(0,exp((-j * ki) + ps));
  new_Matrix_2by2.set(1,0);
  new_Matrix_2by2.set(2,0);
  new_Matrix_2by2.set(3,exp((j * ki) - ps));


  return(new_Matrix_2by2);
}


Tmm::Matrix_2by2 Tmm::get_D(double n1_real,double n1_imag,double n2_real,double n2_imag,double theta_layer1, double theta_layer2)
{
  Tmm::Matrix_2by2 new_Matrix_2by2;
  complex<double> n1_complex (n1_real,n1_imag);
  complex<double> n2_complex (n2_real,n2_imag);
  n1_complex=n1_complex*cos(theta_layer1*M_PI/180);
  n2_complex=n2_complex*cos(theta_layer2*M_PI/180);
  complex<double> r12,r21,t12,t21;

  r12=(n1_complex-n2_complex)/(n1_complex+n2_complex);
  r21=(n2_complex-n1_complex)/(n1_complex+n2_complex);
  t12=(2.0*n1_complex)/(n1_complex+n2_complex);
  t21=(2.0*n2_complex)/(n1_complex+n2_complex);
  new_Matrix_2by2.set(0,1.0/(t12));
  new_Matrix_2by2.set(1,(-r21)/(t12));
  new_Matrix_2by2.set(2,(r12)/(t12));
  new_Matrix_2by2.set(3,(t12*t21-r12*r21)/(t12));

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




vector<double> Tmm::theta_cal(vector<double> n_real , double incident_angle)
{
  vector<double> theta(n_real.size());
  theta[0]=incident_angle;
  for (int k=1; k<n_real.size();k++)
  {
    //theta[k]=asin((n_real[k-1]/n_real[k])*sin(theta[k-1]*M_PI/180))*180/M_PI;
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


    double kr_ratio =0;
    double steps =1;


    // TODO reserve space
    vector<double> n_real;
    vector<double> n_imag;
    vector<double> l_length;
    vector<double> l;
    vector<double> Incoh;


    vector<double> n_real_init;
    vector<double> n_imag_init;
    vector<double> l_length_init;
    vector<double> l_init;
    vector<double> Incoh_init;


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
           if (internal_source_simulation==0 && mod_int->read_type() == "Dipole Source"){
             internal_source_simulation = 1;
             kr_ratio = mod_int->get_kr();
             steps = mod_int->get_steps();
             if (steps == 0)
               steps = 1;
             for (double nm =0; nm<l_init.size();nm++)
               _dipole_loc += l_init[nm];
           }
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
    if (external_source_simulation)
    {
      std::cout<< "Solving for External Source"<<std::endl;
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
      for (int nm = 0; nm<Incoh.size(); ++nm)
        if(Incoh[nm]==1)
          rnd = 5;
      phase_step = 2 * M_PI / rnd;
      double N0 = 0;
      for (int iter = 0; iter <rnd; ++iter )    // loop over added phases
      {
        double randoom = 0;
        double r =_reflectivity;
        //***********************************************************************
        //**********loop added for highly absorbing structures*******************

        Tmm::Matrix_2by2 Unit(1,0,0,1);
        Tmm::Matrix_2by2 DD(1,0,0,1);
        Tmm::Matrix_2by2 MM(0,0,0,0);
        Tmm::Matrix_2by2 TT_load(0,0,0,0);
        Tmm::Matrix_2by2 TT(1,0,0,1);
        for (int k = 0 ;k <= n_real.size()-2 ; ++k)
        {
          DD = get_D(n_real[k],n_imag[k],n_real[k+1],n_imag[k+1],theta[k],theta[k+1]);
          if(Incoh[k]==1 && Incoh[k+1]==0)
          {
            MM = get_M(n_real[k],n_imag[k],l[k],lambda,theta[k],phase_step * iter);
          }
          else
          {
            MM = get_M(n_real[k],n_imag[k],l[k],lambda,theta[k],0);
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
             DD = get_D(n_real[k],-n_imag[k],n_real[k+1],-n_imag[k+1],theta[k],theta[k+1]);
             if(Incoh[k]==1 && Incoh[k+1]==0)
             {
               MM = get_M(n_real[k],n_imag[k],l[k],lambda,theta[k],phase_step * iter);
             }
             else
             {
               MM = get_M(n_real[k],n_imag[k],l[k],lambda,theta[k],0);
             }
             TT_load = DD * T;
             T = MM* TT_load;
             E_I = T * E_N  ;
             E_F.push_back(E_I.get(0)) ;
             E_B.push_back(E_I.get(2));
           }
           //std::cout<<"TMM IS "<<std::endl;
           //T.print();
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


      _Wavelength.push_back(lambda);
      double scale = 1e-6;
      _Reflection.push_back((int)(real(avg_reflection) / scale) * scale);
      _Transmission.push_back((int)(real(avg_transmission) / scale) * scale);
      _Absorption.push_back((int)(avg_Absorption / scale) * scale);

      //****************************************************************************
      //*******************Calculating integral over wavelengths********************
      if (i == 0)
      {

        Electric_Field_integral.resize(n_real.size());
        Electric_Field_integral = Electric_Field_avg;

        intensity_integral.resize(n_real.size());
        intensity_integral = Intensity_avg;

        generation_rate_integral.resize(n_real.size());
        generation_rate_integral = Generation_rate_avg;

      }else
      {
        for (int nm =0; nm< n_real.size();nm++)
        {
          generation_rate_integral[nm] +=  Generation_rate_avg[nm];
          intensity_integral[nm] += Intensity_avg[nm];
          Electric_Field_integral[nm] += Electric_Field_avg[nm];
        }
      }
      if( i == lambda_interp.size()-1)
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
          //solution.add(l_length[nm], generation_rate_integral[nm]);
          _Generation_rate.push_back(generation_rate_integral[nm]);
          _Intensity.push_back(intensity_integral[nm]);
          _External_Source_ElectricField.push_back(Electric_Field_integral[nm]);
        }
      }
    }


    //****************************************************************************
    //*****************************dipole simulation******************************
    //_dipole_loc = 450;
    //internal_source_simulation=1;
    //std::cout<<"dipole loc is "<<_dipole_loc<<std::endl;
    if (internal_source_simulation)
    {
      double dipole_loc = _dipole_loc;
      std::cout<<"Solving for internal source --> dipole loc is "<<dipole_loc<<std::endl;
          double A_P = -sqrt(3/(16 * M_PI));
      double A_N =  sqrt(3/(16 * M_PI));
      double ks = 2*M_PI*n_real[dipole_loc]/lambda;
      for (double kr = 0; kr <= ks*kr_ratio; kr = kr + ks*kr_ratio/steps)
      {
        Tmm::Matrix_2by2 Es(A_P,0,A_N,0);
        Tmm::Matrix_2by2 DD(1,0,0,1);
        Tmm::Matrix_2by2 MM(0,0,0,0);
        Tmm::Matrix_2by2 TT_load(0,0,0,0);
        Tmm::Matrix_2by2 T_RIGHT(1,0,0,1);
        Tmm::Matrix_2by2 T_LEFT(1,0,0,1);
        Tmm::Matrix_2by2 T_TOTAL(1,0,0,1);
        Tmm::Matrix_2by2 AA(0,0,0,1);
          for (double k = n_real.size()-1; k >dipole_loc; --k)
          {
            if (k < n_real.size()-1)
              DD = get_D(n_real[k],n_imag[k],n_real[k+1],n_imag[k+1],0,0);

            MM = get_M(n_real[k],n_imag[k],l[k],lambda,kr,0);
            TT_load = DD * T_RIGHT;
            T_RIGHT = MM* TT_load;
          }
          MM = get_M(n_real[dipole_loc],n_imag[dipole_loc],l[dipole_loc]/2,lambda,kr,0);
          DD = get_D(n_real[dipole_loc],n_imag[dipole_loc],n_real[dipole_loc+1],n_imag[dipole_loc+1],0,0);
          TT_load = DD * T_RIGHT;
          T_RIGHT = MM* TT_load;
          T_LEFT = MM;
          for (double k = dipole_loc-2; k >=0; --k)
          {
            DD = get_D(n_real[k],n_imag[k],n_real[k+1],n_imag[k+1],0,0);
            MM = get_M(n_real[k],n_imag[k],l[k],lambda,kr,0);
            TT_load = DD * T_LEFT;
            T_LEFT = MM* TT_load;
          }
          T_TOTAL = T_LEFT * T_RIGHT;
          AA.set(0, -T_TOTAL.get(0));
          AA.set(2, -T_TOTAL.get(2));
          AA.inv();
          TT_load = T_LEFT * Es;
          TT_load = AA * TT_load;

          /////////////////////////////////////////////////////////////////////////////////
          vector<complex<double>> E_int_f_r(n_real.size()+1);
          vector<complex<double>> E_int_b_r(n_real.size()+1);
          vector<complex<double>> E_int_f_l(n_real.size()+1);
          vector<complex<double>> E_int_b_l(n_real.size()+1);
          vector<complex<double>> E_int_f(n_real.size()+1);
          vector<complex<double>> E_int_b(n_real.size()+1);
          vector<complex<double>> E_int(n_real.size()+1);
          vector<complex<double>> cos_phi;

          E_int_f_r[n_real.size()] = TT_load.get(0);
          E_int_b_r[n_real.size()] = 0;
          E_int_f_l [0] = 0;
          E_int_b_l [0] = TT_load.get(2);
          Tmm::Matrix_2by2 E_N(1,0,0,0);
          Tmm::Matrix_2by2 E_I(0,0,0,0);
          E_N.set(0,E_int_f_r[n_real.size()]);
          E_N.set(2,E_int_b_r[n_real.size()]);
          T_RIGHT.unit_matrix();
          DD.unit_matrix();
          for (double k = n_real.size()-1; k >=dipole_loc; --k)
          {
            if (k < n_real.size()-1)
              DD = get_D(n_real[k],n_imag[k],n_real[k+1],n_imag[k+1],0,0);
            MM = get_M(n_real[k],n_imag[k],l[k],lambda,kr,0);
            TT_load = DD * T_RIGHT;
            T_RIGHT = MM* TT_load;
            E_I = T_RIGHT * E_N;
            E_int_f_r[k] = (E_I.get(0)) ;
            E_int_b_r[k] = (E_I.get(2));
          }
          T_LEFT.unit_matrix();
          E_N.set(2,E_int_f_l[0]);
          E_N.set(0,E_int_b_l[0]);
          DD.unit_matrix();
          for (double k = 0; k <= dipole_loc-1; ++k)
          {
            if (k > 0)
              DD = get_D(n_real[k],n_imag[k],n_real[k-1],n_imag[k-1],0,0);
            MM = get_M(n_real[k],n_imag[k],l[k],lambda,kr,0);
            TT_load = DD * T_LEFT;
            T_LEFT = MM* TT_load;
            E_I = T_LEFT * E_N;
            E_int_f_l[k+1] = (E_I.get(2)) ;
            E_int_b_l[k+1] = (E_I.get(0));
          }

          for (double nm = 0; nm<E_int_f_r.size();nm++)
            if (nm != dipole_loc)
            {
                E_int[nm] = E_int_f_r[nm]+E_int_b_r[nm]+E_int_f_l[nm]+E_int_b_l[nm];
                if (kr == 0)
                {
                _Internal_Source_ElectricField.push_back(real(E_int[nm]));
                _Internal_Intensity.push_back(n_real[nm] * pow(abs(E_int[nm]), 2));
                }
            }
          if (kr == 0)
          {
            _Internal_Source_ElectricField.insert(_Internal_Source_ElectricField.begin()+dipole_loc,((_Internal_Source_ElectricField[dipole_loc]+_Internal_Source_ElectricField[dipole_loc-1])/2));
            _Internal_Intensity.insert(_Internal_Intensity.begin()+dipole_loc,((_Internal_Intensity[dipole_loc]+_Internal_Intensity[dipole_loc-1])/2));
          }
         // for (int yy =0; yy <E_int.size();++yy)
            //std::cout << yy << "   " << E_int[yy] << std::endl;

           // _Output_Front.push_back(n_real[0] * pow(abs(E_int[0]), 2));
            //_Output_back.push_back(n_real[n_real.size()-1] * pow(abs(E_int[E_int.size()-1]), 2));
          vector<double> total_internal_intesity;
          for (double nm = 0; nm<E_int_f_r.size();nm++)
          {
            E_int_f[nm] = E_int_f_l[nm] + E_int_f_r[nm];
            E_int_b[nm] = E_int_b_l[nm] + E_int_b_r[nm];
          }
          E_int_f.erase(E_int_f.begin()+dipole_loc);
          E_int_b.erase(E_int_b.begin()+dipole_loc);
          for (double nm = 0; nm<n_real.size();nm++)
          {
            complex<double> ki ((2*M_PI*n_real[nm]*l[nm])/lambda , (2*M_PI*n_imag[nm]*l[nm])/lambda);
            complex<double> kzi (1,0);
            kzi = ki * sqrt(1-(pow(kr,2)/pow(ki,2)));
            cos_phi.push_back(kzi/ki);
            complex<double> Ni (n_real[nm], n_imag[nm]);
            total_internal_intesity.push_back(real(Ni*cos_phi[nm]*conj(E_int_f[nm]+E_int_b[nm])*(E_int_f[nm]-E_int_b[nm])));
          }
         // std::cout<<"cos phi is "<<cos_phi[0] <<" and " <<cos_phi[total_internal_intesity.size()-1]<<std::endl;
            _Output_Front_angle.push_back(acos(real(cos_phi[0])));
            _Output_Front.push_back(total_internal_intesity[0]);
            _Output_back_angle.push_back(acos(real(cos_phi[total_internal_intesity.size()-1])));
            _Output_back.push_back(total_internal_intesity[total_internal_intesity.size()-1]);
      }
    }// end of loop of dipole simulation


  }// end of loop of wave_length

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
        file << _Wavelength[i] << " " << _Transmission[i] << " "
             << _Reflection[i] << " "
             << _Absorption[i] << "\n";
    }
  }

  file.close();

    string polar_file(outdir + "/" + get_output_filename() + "_polar.dat");
    ofstream polar;
    polar.open(polar_file.c_str());
    if (polar.good())
      polar << '#' << get_type() << " Polar radiation pattern (" << get_name() << ")\n"
            << "# front_angle[deg] " <<  " front_Intensity[a.u.] "
            << "# back_angle[deg] " <<  " back_Intensity[a.u.] " << "\n";
    for (double i = 0; i < _Output_Front.size(); ++i)
      polar << _Output_Front_angle[i] << "   " << _Output_Front[i] << "      " << _Output_back_angle[i]<< "      " <<_Output_back[i] << "\n";

    polar.close();

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

  if (solutions.count(Internal_Intensity))
  {
    solutions[Internal_Intensity][0]= _Internal_Intensity[dof_indices[0]];
  }
  if (solutions.count(Internal_Source_ElectricField))
  {
    solutions[Internal_Source_ElectricField][0]= _Internal_Source_ElectricField[dof_indices[0]];
  }
  if (solutions.count(External_Source_ElectricField))
  {
    solutions[External_Source_ElectricField][0]= _External_Source_ElectricField[dof_indices[0]];
  }
  if (solutions.count(Intensity))
  {
    solutions[Intensity][0]= _Intensity[dof_indices[0]];
  }
  if (solutions.count(GenerationRate))
  {
    solutions[GenerationRate][0] =_Generation_rate[dof_indices[0]];
  }

}

