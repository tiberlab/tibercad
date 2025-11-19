// $Id$

#ifndef _TMM_H_
#define _TMM_H_

#include "tibercad/module/SimulationInterface.h"
#include "matrix2by2.h"
using namespace std;



/*!
 * \brief Implementation of Transfer Matrix Method for electromagnetic fields
 *
 * This class implements standard Transfer Matrix Method (TMM) to calculate the
 * electromagnetic field distribution in a 1D layered structure.
 *
 * Author:
 * Contributors:
 */

class TBDLLOCAL Tmm : public SimulationInterface
{

public:

  //! Destructor
  /*!
   * We do not declare it virtual here, as we will not allow
   * to derive from this class anyway.
   */
  ~Tmm(void);

  //! We need a public static creator function
  static Tmm* create(const ModelOptions& options);



  //! defining a function to return M matrix
  /*!
  * "n_real" is real part of refractive index
  * "n_imag" is imaginary part of refractive index
  * "length" is the lenth of the layer
  * "lambda" is the light's wavelength
  * "theta" is the light's traveling angle(normal incident is equal to '0')
  */
  matrix2by2 get_M(double n_real, double n_imag, double length, double lambda, double kr, double phase);

  matrix2by2 Determinal_Matrix(matrix2by2 MAT);


  std::complex<double> cmlx_sqrt(std::complex<double> in);

  bool cmp_string(const std::string& a, const std::string& b);


  void reset_global_variables(void);

  //! defining a function to return D matrix
  /*!
  * "n1_real" is first layer real part of the refractive index
  * "n1_imag" is first layer imaginary part of the refractive index
  * "n2_real" is second layer real part of the refractive index
  * "n2_imag" is second layer imaginary part of the refractive index
  * "theta_layer1" is the first layer light's traveling angle(normal incident is equal to '0')
  * "theta_layer2" is the second layer light's traveling angle(normal incident is equal to '0')
  */
  matrix2by2 get_D(double n1_real, double n1_imag, double n2_real, double n2_imag, double kr, double lambda, double mode);

  void dipole_source(double& A_P, double& A_N, double Mode, double Oraintation, double lambda, double cos_phi_inter);

  void solving_internal_source(vector<complex<double>>& E_int_f, vector<complex<double>>& E_int_b, const vector<double>& n_real,
    const vector<double>& n_imag, const vector<double>& l, double lambda, double dipole_loc, double kr, double Mode, double A_P, double phase);


  vector<double> linear_interpolation1(vector<double> xData, vector<double> yData, vector<double> x_interp);

protected:


  virtual void plot_globaldata(void);


  //! The initialization
  virtual void do_init(void);


  //! Parse the options from the input file
  virtual void parse_options(void);


  //! Setup the available variables
  virtual void do_setup_solution_variables(void);


  //! Solve the Poisson equation
  virtual void do_solve(void);


  //! Print some useful information
  virtual void do_print_info(void);


  //! We need to create a physical model
  virtual PhysicalModel* create_bulk_model(const ModelOptions& options,
    const Material* mat) const;

  //! We need to create boundary condition model
  virtual PhysicalModel* create_boundary_model(const ModelOptions& options,
    const MaterialBoundary* boundary) const;


  //! We have to provide somehow our solution variables
  virtual void get_solution_secure(const Elem* elem,
    std::map<ID, std::vector<double> >& values,
    const std::vector<Point>& p);

  //! Get a mesh independent solution variable
//  virtual void get_solution_secure(std::map<ID, std::vector<double> >& values);



private:



  //! function to calculate light's propagation angle in all layers
  /*!
   * This function calculate propagating angle of the light using snell's law in all layers
   * "n_real" is a Refractive index vector
   * "incident_angle" is angle of light in the first layer
   * normal incident is equal to incident_angle = 0
   */


  Device* _device;

  SimulationEnvironment* _env;

  //! defining a class to work with 2*2 matrics 



  /*!
   * \brief The known solution variables
   */
  enum Solutions
  {
    GenerationRate,
    Intensity,
    External_Source_Poynting,
    External_Source_ElectricField,
    Internal_Source_ElectricField,
    Internal_Poynting,
    Energy_Loss,
    Internal_Power,
    Internal_Intensity,
    Internal_Absorption,
    Transmission,
    Reflection,
    Absorption,
    Generation_regions,
    Polar,
    AVT
  };

  std::vector<double> _Transmission;
  std::vector<double> _Reflection;
  std::vector<double> _Absorption;
  double _AVT = 0;

  std::vector<vector<double>> _Electric_Field_External;
  std::vector<vector<double>> _Electric_Field_Internal;
  std::vector<vector<double>> _Generation_regions;
  std::vector<vector<double>> _Generation_regions_internal;
  std::vector<std::string> _regions_name;

  std::vector<double> _Intensity;
  std::vector<double> _External_Source_ElectricField;
  std::vector<double> _Wavelength;
  std::vector<double> _Generation_rate;
  std::vector<double> _Poynting_external;
  std::vector<double> _Energy_loss_external;
  std::vector<double> _OutCoupling;

  std::vector<double> _Internal_Wavelength;
  std::vector<double> _Internal_Source_ElectricField;
  std::vector<double> _Internal_Poynting;
  std::vector<double> _Internal_Power;
  std::vector<double> _Internal_Absorption;
  std::vector<double> _Abs;
  std::vector<double> _Internal_Intensity;

  std::vector<double> _angle;
  std::vector<double> _Poynting_front;
  std::vector<double> _Poynting_back;


  /*!
   * \brief Constructor
   *
   * Being private disables further inheritance.
   */
  Tmm(const ModelOptions& options);

  /*!
   * \brief The wavelengths
   */

   /*!
    * \brief The incident angle
    */
  double _incident_angle;
  double _reflectivity;
  double _up_lambda;
  double _down_lambda;
  double _wavelength_steps;
  double _dipole_loc;


  double _coh_mod;

  std::string _polarization;
  std::string _orientation;
  std::vector<double> _polarization_vec;
  std::vector<double> _oraintation_vec;

  double _steps;
  std::vector<double> _ratio;
  std::vector<double> _dipole_power;

  std::vector<double> _wavelength_vector;
  std::vector<double> _dipole_coordinate;




  //! The wavelengths
  std::vector<double> _lambda;


  //! The solar spectrum
  std::vector<double> _spectrum;

  std::vector<double> _eye_value;
  std::vector<double> _eye_wl;


  std::vector<double> _emission_value;
  std::vector<double> _emission_wl;

  std::map<ID, std::vector<double> > _solutions;


  std::vector<std::vector<double>> _green_vector;
  bool _green_vector_solved;









};


#endif // _TMM_H_
