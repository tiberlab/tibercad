// $Id$

#include "ZbRotatedCrystal.h"
#include "Database.h"
#include "Material.h"
#include "Messages.h"
#include "InitFailedException.h"

#include <cmath>

ZbRotatedCrystal::ZbRotatedCrystal(const ModelOptions& options) :
  RotatedCrystal(options)
{
  a_lat = 0;
  set_xyz_mil_direction("x", 1, 0, 0);
  set_xyz_mil_direction("y", 0, 1, 0);
  set_xyz_mil_direction("z", 0, 0, 1);
}



void ZbRotatedCrystal::set_lat_const(const double a) 
{
  a_lat = a;
}

//===============================================================//

void ZbRotatedCrystal::set_xyz_mil_direction(std::string dir, int h, int k, int l)
{
  assert( (dir.compare(std::string("x"))==0) || (dir.compare(std::string("y"))==0) || (dir.compare(std::string("z")) ==0) );

   if (dir.compare("x")==0) {
     x_miller.clear();
 
     x_miller.push_back(h);
     x_miller.push_back(k);
     x_miller.push_back(l);
  };

    if (dir.compare("y")==0) {
      y_miller.clear();

      y_miller.push_back(h);
      y_miller.push_back(k);
      y_miller.push_back(l);
      
  };

  
    if (dir.compare("z")==0) {
      z_miller.clear();

      z_miller.push_back(h);
      z_miller.push_back(k);
      z_miller.push_back(l);
     
  };

}

//===============================================================//

void ZbRotatedCrystal::calculate_lat_consts()
{

  // If 4 Miller indices are given, then we use a hexagonal
  // cell ([111] growth direction)
  
  if (x_miller.size() == 3)
  {
    lat_const_calc[0] = a_lat/sqrt(double( x_miller[0]*x_miller[0] + x_miller[1]*x_miller[1] + x_miller[2]*x_miller[2]));
    lat_const_calc[1] = a_lat/sqrt(double( y_miller[0]*y_miller[0] + y_miller[1]*y_miller[1] + y_miller[2]*y_miller[2]));
    lat_const_calc[2] = a_lat/sqrt(double( z_miller[0]*z_miller[0] + z_miller[1]*z_miller[1] + z_miller[2]*z_miller[2]));
  }
  else // there is no other choice by now
  {
    // calculate hexagonal basis vectors
    double a = a_lat * M_SQRT1_2;
    double c = 1.632993161855452 * a; // sqrt(8/3)

    //Bravais vectors
    Tensor1  Rx(0);
    Rx(1) = a;
    Tensor1  Ry(0);
    Ry(1) = -0.5 * a;
    Ry(2) = 0.866025403784439 * a; // sqrt(3)/2
    Tensor1  Rz(0);
    Rz(3) = c;

    // Reciprocal  basis
    Tensor1 rec_basis1;
    Tensor1 rec_basis2;
    Tensor1 rec_basis3;

    const double volume = Rx * vectorProduct(Ry, Rz);

    rec_basis1 = vectorProduct(Ry, Rz) / volume;
    rec_basis2 = vectorProduct(Rz, Rx) / volume;
    rec_basis3 = vectorProduct(Rx, Ry) / volume;


   lat_const_calc[0] = 1.0 / norm((x_miller[0] - x_miller[2]) * rec_basis1  +
       (x_miller[1] - x_miller[2]) * rec_basis2 +
       x_miller[3]* rec_basis3);

   lat_const_calc[1] = 1.0 / norm((y_miller[0] - y_miller[2]) * rec_basis1 +
       (y_miller[1] - y_miller[2]) * rec_basis2 +
       y_miller[3] * rec_basis3);


   lat_const_calc[2] = 1.0 / norm((z_miller[0] - z_miller[2]) * rec_basis1 +
       (z_miller[1] - z_miller[2]) * rec_basis2 +
       z_miller[3] * rec_basis3);
  }


 //std::cerr << "lat_const_calc  " << lat_const_calc[0] << "   " << lat_const_calc[1] << "    " <<  lat_const_calc[2] << "\n";
 

}

//===============================================================//

void ZbRotatedCrystal::calculate_rot_matrix_miller(std::vector<int> vec_x_mil, std::vector<int> vec_y_mil)
{

  Tensor1 vec_x;
  Tensor1 vec_y;
  Tensor1 vec_z;

  if (vec_x_mil.size() == 3)
  {
    //convert from miller indexes to vectors
    vec_x(1) = a_lat* vec_x_mil[0];
    vec_x(2) = a_lat* vec_x_mil[1];
    vec_x(3) = a_lat* vec_x_mil[2];

    vec_y(1) = a_lat* vec_y_mil[0];
    vec_y(2) = a_lat* vec_y_mil[1];
    vec_y(3) = a_lat* vec_y_mil[2];

    vec_z(1) = a_lat* z_miller[0];
    vec_z(2) = a_lat* z_miller[1];
    vec_z(3) = a_lat* z_miller[2];
    
    bool error = (fabs(vec_x * vec_y) > 1e-9);
    error |= (fabs(vec_x * vec_z) > 1e-9);
    error |= (fabs(vec_y * vec_z) > 1e-9);
    error |= (vec_x * vectorProduct(vec_y, vec_z) < 1e-9);

    if (error)
    {
      Messages m;
      m.error("These crystal directions are inconsistent "
          "(they do not build right-handed orthogonal system):");
      m.indent();
      std::ostringstream os;
      os << "x = [" << x_miller[0] << " " << x_miller[1] << " " << x_miller[2] << "]";
      m.error(os.str());
      os.str("");
      os << "y = [" << y_miller[0] << " " << y_miller[1] << " " << y_miller[2] << "]";
      m.error(os.str());
      os.str("");
      os << "z = [" << z_miller[0] << " " << z_miller[1] << " " << z_miller[2] << "]";
      m.error(os.str());
      os.str("");
      throw InitFailedException("Inconsistent crystal directions");
    }
  }
  else
  {
    double a = a_lat * M_SQRT1_2;
    double c = 1.632993161855452 * a; // sqrt(8/3)

    // Miller basis
    // assign principal directions
    Tensor1 mil1(0);
    mil1(1) = a;
    Tensor1 mil2(0);
    mil2(1) = -0.5 * a;
    mil2(2) = 0.866025403784439 * a; // sqrt(3)/2
    Tensor1 mil3(0);
    mil3(3) = c;

    //convert from miller indexes to vectors
    vec_x = (x_miller[0] - x_miller[2]) * mil1 +
        (x_miller[1] - x_miller[2]) * mil2 +
        x_miller[3] * mil3;
    vec_y = (y_miller[0] - y_miller[2]) * mil1 +
        (y_miller[1] - y_miller[2]) * mil2 +
        y_miller[3] * mil3;
  }

  calculate_rot_matrix(vec_x, vec_y);
}


//====================================================//
void ZbRotatedCrystal::read_database ( )
{
  const Database& db = get_database();
  db.set_section("lattice");

  //a_lat = db.get("a", 0.543095, true);
  a_lat = db.get("a", 0.543095);
  
   
}

//=====================================================//

void ZbRotatedCrystal::do_init(void)
{

  ModelOptions & options = get_options();
  a_lat = options.get_option("a", a_lat);
  assert(a_lat > 0);

  options.get_option("x-growth-direction", x_miller);
  options.get_option("y-growth-direction", y_miller);
  options.get_option("z-growth-direction", z_miller);

  const std::string& mat = get_material()->get_name();
  int miller_size = x_miller.size();
  if (y_miller.size() != miller_size ||
      z_miller.size() != miller_size)
    throw InitFailedException(mat + ": zincblende growth directions are wrong."
        " (Different number of Miller indices for different directions)");

  if (miller_size == 4)
    Messages::info("Using hexagonal lattice for material " + mat);
  else if (miller_size != 3)
    throw InitFailedException(mat + ": zincblende growth directions are wrong."
        " Can only take 3 or 4 Miller indices.");

  calculate_lat_consts();

  calculate_rot_matrix_miller(x_miller, y_miller);

  calculate_euler_angles();

}

//-----------------------------------------------------//

void  ZbRotatedCrystal::copy_from (const PhysicalModel *rhs)
{

  const ZbRotatedCrystal* temp = dynamic_cast<const ZbRotatedCrystal*> (rhs);
  a_lat = temp->a_lat;

  x_miller = temp->x_miller;
  y_miller = temp->y_miller;
  z_miller = temp->z_miller;


  calculate_lat_consts();

  calculate_rot_matrix_miller(x_miller, y_miller);
  
}

//-----------------------------------------------------//

void  ZbRotatedCrystal::do_init_alloy (const PhysicalModel *comp_A, 
                                       const PhysicalModel *comp_B, double xa) 
{
  
  const ZbRotatedCrystal* mod_A = dynamic_cast<const ZbRotatedCrystal*> (comp_A);
  const ZbRotatedCrystal* mod_B = dynamic_cast<const ZbRotatedCrystal*> (comp_B);

  a_lat = alloy(mod_A->a_lat, mod_B->a_lat, xa);


  x_miller = mod_A->x_miller;
  y_miller = mod_A->y_miller;
  z_miller = mod_A->z_miller;


  calculate_lat_consts();

  calculate_rot_matrix_miller(x_miller, y_miller);
}

//==================================================================//
PhysicalModel* ZbRotatedCrystal::create_new(void) const
{
  return new  ZbRotatedCrystal(get_options());

}

//==================================================================//

