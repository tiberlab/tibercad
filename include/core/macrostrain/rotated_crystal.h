#include "tensor.h"
#include <cmath>
#include <vector>
using namespace std;


class rotated_crystal
{
 public:
  rotated_crystal();
  rotated_crystal(string type, double a);
  rotated_crystal(string type, double a, double c);
  Tensor2Gen RotMatrix; //Rotation matrix  CALC_system_vector_{i} = RotMatrix_{ij}*cryst_system_vector_{j}

                        /*In matrix notation:  
			  vector transformation:    vec_calc = R * vec_crystal
                          2D tensor transformation: Matr_calc = R * Mat_crystal * (R^T)
			 */

  //If one needs to trasform from calculation to crystal, then R^T has to be used
  //------------------------------------------------------------------
  //set and get crystal tyme
  void set_cryst_type( string type_crystal ); 
  string get_cryst_type();
  //-------------------------------------------------------------------
  //sets lattice constant parameters
  void set_lat_const(const double a);
  void set_lat_const(const double a, const double c);
    
  //------------------------------------------------------------------
  //define principal directions by use of Miller indexes
  void set_xyz_mil_direction(string dir, int h, int k, int l);
  void set_xyz_mil_direction(string dir, int h, int k, int i, int l);
  //-------------------------------------------------------------------
 
  void calculate_lat_consts();// calculates Bravais vectors, rotation matrix and lattice constants in calculation system 

  void get_lat_const(double lat_cont[3]); // returns a_lat_calc
  double  lat_const_calc[3];
 
  Tensor2Sym get_eps0(double lat_cont_substrate[3]);

  Tensor2Sym get_const_eps0(double lat_cont_substrate[3], Tensor2Sym& eps0_var_log);

  Tensor2Sym get_var_eps0(std::string  var_name);

 
  void calculate_rot_matrix(const Tensor1& vec_x,const Tensor1& vec_y);

  void calculate_rot_matrix(std::vector<int> vec_x_mil, std::vector<int> vec_y_mil);


  
 private:

  string  crystal_type;  //"cub" - cubic, "hex" - hexagonal
  double  a_lat   ; // lattice constants of crystal
  double  c_lat   ; 
  
  //---------------------------------------------------------
  //Miller indexes of the growth direction

 
  std::vector<int> x_miller; 
  std::vector<int> y_miller;
  std::vector<int> z_miller;
  
 

 
};

