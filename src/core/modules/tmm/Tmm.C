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



Tmm::Tmm(const ModelOptions& options) :
  SimulationInterface(options),
  _incident_angle({0.0})
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
  }
}




Tmm::Matrix_2by2 Tmm::get_M(double n_real,double n_imag,double lenght,double lambda, double theta, double phase)
{
  Tmm::Matrix_2by2 new_Matrix_2by2;
  complex<double> bi ((2*M_PI*n_imag*lenght)/lambda , (2*M_PI*n_real*lenght)/lambda);
  bi=bi* cos(theta*M_PI/180);

  complex<double> ps (0,phase);

  new_Matrix_2by2.set(0,exp(bi + ps));
  new_Matrix_2by2.set(1,0);
  new_Matrix_2by2.set(2,0);
  new_Matrix_2by2.set(3,exp(-bi - ps));


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
  system.add_variable("E", libMeshEnums::CONSTANT, MONOMIAL, &get_region_ids());

 
  system.init();
}


void
Tmm::parse_options(void)
{
  // read wavelengths from input

  get_option("incident_angle", _incident_angle);
  if (_incident_angle.empty())
  {
    Messages::warning("You did not provide any incident_angle for TMM.");
  }


  get_option("reflectivity", _reflectivity);
  if (_reflectivity.empty())
  {
    Messages::warning("You did not provide any reflectivity for TMM.");
  }

  get_option("up_lambda", _up_lambda);
  if (_up_lambda.empty())
  {
    Messages::warning("You did not provide any up_lambda for TMM.");
  }
  get_option("down_lambda", _down_lambda);
  if (_down_lambda.empty())
  {
    Messages::warning("You did not provide any down_lambda for TMM.");
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
  declare_solution(Intensity, REAL, CELL, "W/m^3/nm");
  //declare_solution(HField, VECTOR, CELL, "A/cm");
  //declare_solution(Displacement, VECTOR, CELL, "C/cm^2");
  //declare_solution(Displacement, VECTOR, CELL, "C/cm^2");

  // we can define aliases but the same name cannot refer to
  // several IDs
  //add_alias("ElectricField", EField);
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
  vector<unsigned int> dof_indices;
  vector<double> area;
  vector<double> gen;
  vector<double> lambda_interp;
  vector<double> sun_interp;
  for (double i = _down_lambda[0]; i<= _up_lambda[0]; i += 1)
  {
    lambda_interp.push_back(i);
  }
  sun_interp = Tmm::linear_interpolation1(_lambda,_spectrum,lambda_interp);


  for (unsigned int i = 0; i < lambda_interp.size(); ++i)   //loop over wavelength
  {

  for (unsigned int j = 0; j < _incident_angle.size(); ++j)  //loop over incident angle
  {

    std::cout<<"solving condition :  " << std::endl;
    std::cout<<"Lambda is  " << lambda_interp[i] <<  std::endl;
    std::cout<<"Radiation intensity is  " << sun_interp[i] <<  std::endl;
   // std::cout<<"Incident angle is  " << _incident_angle[j] <<  std::endl;

    double lambda = lambda_interp[i];



    const unsigned int uvar = system.variable_number("E");


    // TODO this will not work if the 1D mesh is distributed. In that case, MPI calls could be used
    // to gather pieces from all processes
    MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();

    double E = 1.0;
    double incoming_angle=_incident_angle[j];

    double c0 = 2.998e8 * 1e9;              // Speed of light [nm/s]
    double e0 = 8.85e-12;                   // Vacuum permittivity  [F/m]
    double radiation = sun_interp[i] / 1e3;                   //[W/(m^2*nm)]
    double Esun2 = 2 * radiation / (c0 * 1e-9 * e0); // E^2 Electric Field [V^2 / (m^2 * nm)]
    double w = 2 * M_PI * c0 / lambda;                               // Frequency [1/s]
    double plank_const = 1.055e-34;

    // TODO reserve space
    vector<double> n_real;
    vector<double> n_imag;
    vector<double> l_length;
    vector<double> l;


    vector<double> n_real_init;
    vector<double> n_imag_init;
    vector<double> l_length_init;
    vector<double> l_init;


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
      
      // getting refractive index and length of layers from model
      libMesh::Complex nk = mod.get_refractive_index(lambda);
      // n_real.push_back(real(nk));
      // n_imag.push_back(imag(nk));

      n_real_init.push_back(real(nk));
      n_imag_init.push_back(imag(nk));


      //n_imag.push_back(sqrt((abs(mod.get_permittivity(lambda))-real(mod.get_permittivity(lambda)))/2));
      //n_real.push_back(sqrt((abs(mod.get_permittivity(lambda))+real(mod.get_permittivity(lambda)))/2));
      //l.push_back(elem->volume());
      l_init.push_back(elem->volume());
     // std::cout<<"elem n_sides is :" << elem->n_sides() << std::endl;

      // the sides

       for (unsigned int s = 0; s < elem->n_sides(); s++)
       {
         TmmBoundaryModel* mod_int =
           get_interface_model<TmmBoundaryModel>(elem, s);

        // std::cout<<"elements  "<<s<<"   "<<elem->n_sides()<<std::endl;
         if (mod_int != NULL)
         {
           if (mod_int->read_type() == "Mirror"){
             mod_int->Calculate_M_Matrix();
             BC_check_init.push_back(1);
             Boundry_condition.set(0,mod_int->get_element(0));
             Boundry_condition.set(1,mod_int->get_element(1));
             Boundry_condition.set(2,mod_int->get_element(2));
             Boundry_condition.set(3,mod_int->get_element(3));
           }

           if (mod_int->read_type() == "Incident Wave"){
             BC_check_init.push_back(0);
             if ( s % 2 == 0)
             {
               direction = "top to down propagation";
               //std::cout<<direction<<std::endl;
             }
             else
             {
               direction = "down to top propagation";
               //std::cout<<direction<<std::endl;
             }
           }

         }
           if (mod_int == NULL)
           {
             BC_check_init.push_back(0);
           }
       }
    }


    if (direction == "top to down propagation")
    {
      for (int uu=1; uu<BC_check_init.size();uu += 2)
        BC_check.push_back(BC_check_init[uu]);
      for (int uu=0; uu<n_real_init.size();++uu)
      {
        n_real.push_back(n_real_init[uu]);
        n_imag.push_back(n_imag_init[uu]);
        l.push_back(l_init[uu]);
      }

    }else if (direction == "down to top propagation          ")
    {
      for (int uu=BC_check_init.size()-2; uu>=0;uu -= 2)
        BC_check.push_back(BC_check_init[uu]);
      for (int uu=n_real_init.size()-1; uu>=0;--uu)
      {
        n_real.push_back(n_real_init[uu]);
        n_imag.push_back(n_imag_init[uu]);
        l.push_back(l_init[uu]);
      }
    }


    //********************snell's law********************
    vector<double> theta(n_real.size());
    theta=Tmm::theta_cal(n_real,incoming_angle);



    //****************************************************
    vector<double> phase(n_real.size());
    srand( time( NULL ) );  // Initialize random seed
    for(int uu =0; uu<phase.size(); ++uu)
    {
      phase[uu] =  2*M_PI*(double) rand()/RAND_MAX;      // Pseudo-random number between 0 and 1 * 2*Pi

    }



    //*****************************************************
    //************** defining Vectors**********************
    
    Tmm::Matrix_2by2 D(1,0,0,1);
    Tmm::Matrix_2by2 M(0,0,0,0);
    Tmm::Matrix_2by2 T_load(0,0,0,0);
    Tmm::Matrix_2by2 T(1,0,0,1);
    
    double r =_reflectivity[0];
    Tmm::Matrix_2by2 E_N(1,0,r,0);
    Tmm::Matrix_2by2 E_I(0,0,0,0);


    vector<complex<double>> E_F(n_real.size());
    vector<complex<double>> E_B(n_real.size());
    
   // E_N.print();

    E_F[n_real.size()-1]=E_N.get(0);
    E_B[n_real.size()-1]=E_N.get(2);


    //******************************************************
    //******main loop over layer, calculating matrixs*******
    for (int k=n_real.size()-1 ; k>=0 ; --k){
      if (BC_check[k] == 0) {
        if (k<n_real.size()-1){
         D = get_D(n_real[k],n_imag[k],n_real[k+1],n_imag[k+1],theta[k],theta[k+1]);
         //std::cout<<"BC_check[k] == 0"<< std::endl;
        }

      }else {
        D.set(0,Boundry_condition.get(0));
        D.set(1,Boundry_condition.get(1));
        D.set(2,Boundry_condition.get(2));
        D.set(3,Boundry_condition.get(3));
        //std::cout<<"BC_check[k] == 1"<< std::endl;
      }
      //std::cout<<"D matrix is :" <<k<<"   "<< std::endl;
      //D.print();
     // M = get_M(n_real[k],n_imag[k],l[k],lambda,theta[k],phase[k]);
      M = get_M(n_real[k],n_imag[k],l[k],lambda,theta[k],0);
      //std::cout<<"M matrix is :" <<k<<"   "<< std::endl;
      //M.print();
      T_load = D * M;
      T = T * T_load;

      if (k<n_real.size()-1){
      E_I = T * E_N ;

      //std::cout<<"E_I is "<<k <<std::endl;
      //E_I.print();
      E_F[k]+= E_I.get(0);
      E_B[k]+= E_I.get(2);
      }
     // std::cout<<"E_F is "<<k<<"    "<<abs(E_F[k]) <<std::endl;
    }
    cout<<"TMM is :"<<endl;
    T.print();

    //***********************************************************************
    //**********reflection and transmission calculation**********************

    complex<double> Reflection,Transmission;
    Reflection = pow(abs(T.get(2)/T.get(0)),2);
    complex<double> nc_first (n_real[0],n_imag[0]);
    complex<double> nc_last (n_real[n_real.size()-1],n_imag[n_imag.size()-1]);
    complex<double> ratio_complex;
    ratio_complex = ((nc_last)*cos(theta[theta.size()-1]*M_PI/180))/(nc_first*cos(theta[0]*M_PI/180));
    Transmission = ratio_complex*pow(abs(1.0/T.get(0)),2);

    //***************************************************************************
    //**************printing tansmission and reflection**************************
    cout<<"trasmision is :"<< Transmission << "reflection is :"<<Reflection<<endl;
    cout<<"Sum is :"<< Transmission+Reflection<<endl;

    //****************************************************************************
    //***************normalizing electric field matrix******************************

    vector<complex<double>> E_F_NORM(n_real.size());
    vector<complex<double>> E_B_NORM(n_real.size());

    for(double nm=0; nm < n_real.size() ; ++nm)
    {
      E_F_NORM[nm] = E_F[nm]/E_F[0];
      E_B_NORM[nm] = E_B[nm]/E_F[0];
    } 







    vector<complex<double>> Etot(n_real.size()); // Electric Field, Magnetic Field
    vector<complex<double>> Intensity;// Intensity, forward, backward
    vector<complex<double>> Generation_rate;        // Generation Rate
    vector<double> Generation_rate_real;        // Generation Rate real



    for (int nm=0 ; nm < n_real.size() ; ++nm)
    {
      Etot[nm] = E_F_NORM[nm] + E_B_NORM[nm];
      Intensity.push_back(0.5 * c0 * 1e-9 * e0 * n_real[nm] * Esun2* pow(abs(Etot[nm]), 2)); // Intensity [W/(m^2 * nm)]
      Generation_rate.push_back(1 / (plank_const* w) * (4 * M_PI * n_imag[nm] * 1e7/(lambda)) * real(Intensity[nm])/ 1e4 );        // Generation rate [ cm^-3 s^-1 nm^-1]
      Generation_rate_real.push_back(real(Generation_rate[nm]));
    }




    // taking integral over the bandwidths

    if (i == 0)
    {
      area.resize(n_real.size());
      area = Generation_rate_real;
    }

    if (i != 0 && i != lambda_interp.size()-1)
    {
      for (int nm =0; nm< n_real.size();nm++)
         area[nm] +=2 * Generation_rate_real[nm];
    }
    if (i == lambda_interp.size()-1)
    {
      for (int nm =0; nm< n_real.size();nm++)
         area[nm] += Generation_rate_real[nm];
      for (int nm =0; nm< n_real.size();nm++)
         area[nm] *= (lambda_interp[lambda_interp.size()-1] - lambda_interp[0]) / (2 * lambda_interp.size()-1);

    }


    if( i == lambda_interp.size()-1)
    for(double nm=0; nm < l_length.size() ; ++nm){
      //solution.add(l_length[nm], real(Generation_rate[nm]));
      //solution.add(l_length[nm], real(Intensity[nm]) );
      solution.add(l_length[nm], area[nm]);
    //  std::cout<<"area is " << nm << "   " << area[nm] <<std::endl;
    }






    //*******************************************************************************
    //************************printing electric field********************************

    /////////////////////////////////////////////////////////
   // for(int nm=0; nm<_lambda.size();++nm){
      //std::cout<<"sun radiation is " << nm << "   " << _lambda[nm] << "   "<< _spectrum[nm] <<std::endl;
   // }





    //////////////////////////////////////////////////////////////////



  }// end of loop of angle


  }// end of loop of wave_length



  solution.close();
  system.update();
  //solution.print_matlab("sol.m");

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
Tmm::get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p)
{
  unsigned int np = p.size();

  TiberLinearSystem& system = get_equation_system<TiberLinearSystem>();

  const NumericVector<libMesh::Number>& solution = system.get_solution_vector();

  const unsigned int dim = get_mesh().mesh_dimension();

  const DofMap& dof_map = system.get_dof_map();

  const unsigned int u_var = system.variable_number("E");

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


  for (unsigned int n = 0; n < np; n++)
  {
    //double efield  = 0.0;

    //if (values.count(EField))
    //  values[EField][n] = efield;

  }


  if (values.count(GenerationRate))
  {
    values[GenerationRate][0] = solution(dof_indices[0]);
  //  values[EField][1] = field(1) / np;
  //  values[EField][2] = field(2) / np;
  }

}

