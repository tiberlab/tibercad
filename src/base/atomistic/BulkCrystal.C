#include "BulkCrystal.h"
#include "CrystalDefs.h"
#include "Database.h"
#include "Alloy.h"
#include "RuntimeException.h"

#include <boost/algorithm/string/replace.hpp>

#include <fstream>
#include <numeric>

BulkCrystal*
BulkCrystal::create(const Material* mat, const ModelOptions& options)
{
  BulkCrystal* bc = nullptr;
  if ((mat != nullptr) && (mat->get_structure() != "none"))
    bc = new BulkCrystal(mat, options);

  return bc;
}

BulkCrystal::BulkCrystal(const Material* mat, const ModelOptions& options)
{
  _options = ModelOptions(mat->get_options()); 
  _options += options;
  _mat = mat;
}

void
BulkCrystal::init(void)
{
 
  //Careful with the order of the calls, it's important
  read_database();

  set_cell_vectors();

  _basis = _lattice_basis;
  for (auto it = _basis.begin(); it != _basis.end(); ++it)
  {
    it->set_position(_prim_vec * it->get_ttype_position());
  }

  //! Get rotation and apply to basis and primitive vectors
  build_rotation();
  set_ttype_lattice_vectors(_rotated_prim_vec);
  _atoms = _rotated_basis;
  
  // calculate the Euler angles for the rotation matrix
  calculate_euler_angles();

  Messages m;
  m.indent();

  // write some info on the atomic structure
  /*
  std::ostringstream os;
  os << "Basis atom | label\n";
  os << "==================\n";

  for (auto&& a : _basis)
  {
    os << std::left << std::setw(11) << a.get_specie() << "|"
       << std::right << std::setw(6) << static_cast<unsigned int>(a.get_label())
       << "\n";
  }
  m.info(os.str());
  m.newline();
  */
   

  build_bond_map();
  refresh();
  
}


bool
BulkCrystal::extract_crystal_direction(const std::string& dir, Tensor1& vec) const
{
  bool valid = false;

  vec = 0;

  std::string opt_name = dir + "-growth-direction";

  std::string miller_str = _options.get_option(opt_name, "");

  if (!miller_str.empty())
  {
    // the Miller (or Miller-Bravais) indices, can be for directions or planes,
    // depending on type of parentheses
    std::vector<int> miller;
    Utils::extract_vector(miller_str, miller);

    if ((miller.size() > 2) && (miller.size() < 5))
    {
      Tensor1 mv;
      mv(1) = miller[0];
      mv(2) = miller[1];
      mv(3) = miller.back();

      if (miller_str.front() == '[')
      {
        // identify Miller ind. as crystal directions

        if (miller.size() == 4)
        {
          mv(1) = (2 * miller[0] + miller[1]);
          mv(2) = (2 * miller[1] + miller[0]);
        }

        vec = _conv_vec * mv;
      }
      else
      {
        // identify Miller ind. as crystal planes
        vec = _reciprocal_lattice[0] * mv(1) +
              _reciprocal_lattice[1] * mv(2) +
              _reciprocal_lattice[2] * mv(3);
      }

      valid = true;
    }
    else
    {
      std::ostringstream os;
      os << "Wrong Miller indices: " << miller_str;
      throw InitFailedException(os.str()); 
    }
  }

  return(valid);
}

  
void
BulkCrystal::get_orthogonal_vector(const Tensor1& dir, Tensor1& ortho) const
{
  ortho = Tensor1(1.0);

  if (abs(dir(1)) > 1e-3)
  {
    ortho(1) = -(dir(2) + dir(3)) / dir(1);
    return;
  }

  if (abs(dir(2)) > 1e-3)
  {
    ortho(2) = -(dir(1) + dir(3)) / dir(2);
    return;
  }

  if (abs(dir(3)) > 1e-3)
  {
    ortho(3) = -(dir(1) + dir(2)) / dir(3);
    return;
  }
  
}

void
BulkCrystal::build_rotation(void)
{
  std::ostringstream os;

  // reset to unit element
  _rotation = Tensor2Gen(1);

  // the crystal directions along x, y, z
  Tensor1 vec_x;
  Tensor1 vec_y;
  Tensor1 vec_z;

  bool has_x;
  bool has_y;
  bool has_z = false;

  bool has_euler = false;

  std::vector<double> euler_angles;
  _options.get_option("euler_angles", euler_angles);
  if (euler_angles.size() == 3)
  {
    has_euler = true;

    double a = M_PI * euler_angles[0] / 180.0;
    double b = M_PI * euler_angles[1] / 180.0;
    double c = M_PI * euler_angles[2] / 180.0;

    // calculate rotation matrix
    double ca = cos(a);
    double cb = cos(b);
    double cc = cos(c);
    double sa = sin(a);
    double sb = sin(b);
    double sc = sin(c);

    _rotation(1,1) = ca*cb*cc - sa*sc;
    _rotation(1,2) = sa*cb*cc + ca*sc;
    _rotation(1,3) = sb*cc;
    _rotation(2,1) = -ca*cb*sc - sa*cc;
    _rotation(2,2) = ca*cc - sa*cb*sc;
    _rotation(2,3) = -sb*sc;
    _rotation(3,1) = -ca*sb;
    _rotation(3,2) = -sa*sb;
    _rotation(3,3) = cb;

    // get rotated x and y
    has_y = has_x = true;
    vec_x(1) = _rotation(1,1);
    vec_x(2) = _rotation(1,2);
    vec_x(3) = _rotation(1,3);
    vec_y(1) = _rotation(2,1);
    vec_y(2) = _rotation(2,2);
    vec_y(3) = _rotation(2,3);
  }
  else
  {
    if (euler_angles.size() > 0)
      throw InitFailedException("All 3 Euler angles have to be specified.");

    has_x = extract_crystal_direction("x", vec_x);
    has_y = extract_crystal_direction("y", vec_y);
    has_z = extract_crystal_direction("z", vec_z);
  }
  


  // if all are given, make some sanity checks
  if (has_x && has_y && has_z)
  {
    double v = vec_x * vectorProduct(vec_y, vec_z);
    if (v < 1e-12)
    {
      throw InitFailedException("The given crystal directions are inconsistent "
              "(they do not build right-handed system)");
    }
  }

  if (has_x)
  {
    if (has_y && !has_z)
      vec_z = vectorProduct(vec_x, vec_y);

    if (!has_y && has_z)
      vec_y = vectorProduct(vec_z, vec_x);

    if (!has_y && !has_z)
    {
      get_orthogonal_vector(vec_x, vec_y);
      vec_z = vectorProduct(vec_x, vec_y);
    }
  }
  else if (has_y)
  {
    if (has_z)
    {
      vec_x = vectorProduct(vec_y, vec_z);
    }
    else
    {
      get_orthogonal_vector(vec_y, vec_x);
      vec_z = vectorProduct(vec_x, vec_y);
    }
  }
  else if (has_z)
  {
    get_orthogonal_vector(vec_z, vec_y);
    vec_x = vectorProduct(vec_y, vec_z);
  }

  {
    double a = vec_x * vec_y;
    double b = vec_x * vec_z;
    double c = vec_y * vec_z;
    if (((a * a) + (b * b) + (c * c)) > 1e-9)
    {
      throw InitFailedException("The given crystal directions are inconsistent "
                                "(vectors are non-orthogonal)");
    }
  }

  if (!has_euler)
  {
    for (int i = 1; i <= 3; i++)
    {
      _rotation(1, i) = vec_x(i) / norm(vec_x);
      _rotation(2, i) = vec_y(i) / norm(vec_y);
      _rotation(3, i) = vec_z(i) / norm(vec_z);
    }
  }


  _rotated_prim_vec = _rotation * _prim_vec;
  _rotated_conv_vec = _rotation * _conv_vec;

  //! Keep a copy of rotated crystal basis
  _rotated_basis = _basis;
  for (std::vector<Atom>::iterator it = _rotated_basis.begin();
       it != _rotated_basis.end(); ++it)
  {
    it->set_position(_rotation * it->get_ttype_position());
  }

  // recalculate Miller indices and make sanity check
  // we also set the options to the calculated values
  std::vector<int> mil_x;
  get_miller_indices(vec_x, mil_x);
  std::vector<int> mil_y;
  get_miller_indices(vec_y, mil_y);
  std::vector<int> mil_z;
  get_miller_indices(vec_z, mil_z);

  bool miller_bravais = false;
  if (has_x && !has_euler)
  {
    std::string miller_str = _options.get_option("x-growth-direction", "");
    std::vector<int> miller;
    Utils::extract_vector(miller_str, miller);

    miller_bravais |= (miller.size() == 4);

    if (miller_str.front() == '(')
    {

      if ((miller[0] != mil_x[0]) ||
          (miller[1] != mil_x[1]) ||
          (miller.back() != mil_x[2]))
      {
        std::ostringstream os;
        os << "Calculated Miller indices are not the same as the ones provided: "
           << miller_str << " != " << "(" << mil_x[0] << mil_x[1];
        (miller.size() == 4) ? os << -(mil_x[0] + mil_x[1]) : os << "";
        os << mil_x.back() << ")";

        throw InitFailedException(os.str());
      }
    }
  }
  if (has_y && !has_euler)
  {
    std::string miller_str = _options.get_option("y-growth-direction", "");
    std::vector<int> miller;
    Utils::extract_vector(miller_str, miller);

    miller_bravais |= (miller.size() == 4);

    if (miller_str.front() == '(')
    {

      if ((miller[0] != mil_y[0]) ||
          (miller[1] != mil_y[1]) ||
          (miller.back() != mil_y[2]))
      {
        std::ostringstream os;
        os << "Calculated Miller indices are not the same as the ones provided: "
           << miller_str << " != " << "(" << mil_y[0] << mil_y[1];
        (miller.size() == 4) ? os << -(mil_y[0] + mil_y[1]) : os << "";
        os << mil_y.back() << ")";

        throw InitFailedException(os.str());
      }
    }
  }
  if (has_z && !has_euler)
  {
    std::string miller_str = _options.get_option("z-growth-direction", "");
    std::vector<int> miller;
    Utils::extract_vector(miller_str, miller);

    miller_bravais |= (miller.size() == 4);

    if (miller_str.front() == '(')
    {

      if ((miller[0] != mil_z[0]) ||
          (miller[1] != mil_z[1]) ||
          (miller.back() != mil_z[2]))
      {
        std::ostringstream os;
        os << "Calculated Miller indices are not the same as the ones provided: "
           << miller_str << " != " << "(" << mil_z[0] << mil_z[1];
        (miller.size() == 4) ? os << -(mil_z[0] + mil_z[1]) : os << "";
        os << mil_z.back() << ")";

        throw InitFailedException(os.str());
      }
    }
  }

  // now we reconstruct the vectors
  vec_x = mil_x[0] * _reciprocal_lattice[0] + 
          mil_x[1] * _reciprocal_lattice[1] +
          mil_x[2] * _reciprocal_lattice[2];
  
  vec_y = mil_y[0] * _reciprocal_lattice[0] + 
          mil_y[1] * _reciprocal_lattice[1] +
          mil_y[2] * _reciprocal_lattice[2];

  vec_z = mil_z[0] * _reciprocal_lattice[0] + 
          mil_z[1] * _reciprocal_lattice[1] +
          mil_z[2] * _reciprocal_lattice[2];


  // calculate the lattice constants along calculation system axes
  // However, I'm not sure they really serve something
  _ortho_lattice_constants[0] = 1.0 / norm(vec_x);
  _ortho_lattice_constants[1] = 1.0 / norm(vec_y);
  _ortho_lattice_constants[2] = 1.0 / norm(vec_z);


  // set the module options to the calculated directions
  // this is mainly for printing info afterwards
  if (miller_bravais)
  {
    mil_x.insert(mil_x.end()-1, -(mil_x[0] + mil_x[1]));
    mil_y.insert(mil_y.end()-1, -(mil_y[0] + mil_y[1]));
    mil_z.insert(mil_z.end()-1, -(mil_z[0] + mil_z[1]));
  }
  _options.set_option("x-growth-direction", mil_x);
  _options.set_option("y-growth-direction", mil_y);
  _options.set_option("z-growth-direction", mil_z);
}



void
BulkCrystal::calculate_strain(Tensor2Gen &strain,
                 const Tensor2Gen &reference,
                 const Tensor2Gen &other_cell) const
{
  strain = inv(other_cell) * reference - Tensor2Gen(1.0);
}



void
BulkCrystal::find_least_common_lattice(const BulkCrystal &substrate,
                                       libMesh::RealTensor &trafo,
                                       libMesh::RealTensor &residual_strain,
                                       std::vector<unsigned int> &super_ref,
                                       std::vector<unsigned int> &super_mat,
                                       double max_strain) const
{

  double minResidualStrain = std::numeric_limits<double>::max();

  Tensor2Gen inv_prim = inv(_rotated_prim_vec);

  // Loop over possible supercell multiples (up to some reasonable limit)
  int maxSupercell = 1; // Let's limit to multiples up to 3 for simplicity

  Tensor2Gen diag(1);
  Tensor2Gen superA;
  Tensor2Gen superB;
  Tensor2Gen strain(0.0);
  Tensor2Gen deformation(1.0);

  double residualStrain;

  int i = 1, j, k;
  for ( ; i <= maxSupercell; ++i)
  {
    diag(1, 1) = i;
    for (j = 1; j <= maxSupercell; ++j)
    {
      diag(2, 2) = j;
      for (k = 1; k <= maxSupercell; ++k)
      {
        diag(3, 3) = k;

        // Form the supercells for both crystals
        superA = diag * substrate.get_rotated_prim_vec();

        // Try to find a linear transformation that maps our cell to superA
        deformation = inv_prim * superA;

        // extend to least integers
        //Utils::scale_to_int(deformation, 0.02);

        superB = _rotated_prim_vec * deformation;

        // Compute the residual strain tensor
        calculate_strain(strain, superA, superB);
        residualStrain = trace(strain);

        // If residual strain is within tolerance, update the best deformation
        if (abs(residualStrain) < minResidualStrain)
        {
          minResidualStrain = abs(residualStrain);

          if (abs(residualStrain) < max_strain)
            break;
        }
      }
      if (abs(residualStrain) < max_strain)
        break;
    }
    if (abs(residualStrain) < max_strain)
      break;
  }
/*
  std::cerr << i << " " << j << " " << k << "\n";
  std::cerr << superA;
  std::cerr << superB;
  std::cerr << deformation;
  std::cerr << residualStrain << "\n";
*/
}



bool
BulkCrystal::get_lattice_matching_strain(const BulkCrystal& substrate,
                                  libMesh::RealTensor& strain) const
{
  bool compatible = (get_lattice_type() == substrate.get_lattice_type());

  // TODO check also zb (111) on wz (0001) ore vice versa

  libMesh::RealTensor trafo;
  libMesh::RealTensor residual_strain;
  std::vector<unsigned int> super_ref;
  std::vector<unsigned int> super_mat;
  find_least_common_lattice(substrate, trafo, residual_strain,
  super_ref, super_mat, 0.001);

  if (compatible)
  {
    // this one is strained to bulk crystal, and strain tensor rotated
    // to calculation system.

    // cartesian basis vectors in terms of conventional cell are
    // (x y z) = inv(_conv_vec)*_conv_vec
    // so we calculate strained (x' y' z') by using
    // A = inv(_conv_vec) of the strained layer and
    // for _conv_vec the one of the substrate, i.e. we strain the 
    // conventinal vectors and calculate the according strain in 
    // cartesian coordinates:
    // (x', y', z')' = (1 + eps)*(x, y, z)'
    // => eps = (x', y', z') - 1, because (x, y, z) = 1

    Tensor2Gen A(inv(_conv_vec));

    Tensor2Gen eps = A * substrate.get_conv_vec() - Tensor2Gen(1.0);

    // symmetrize to eliminate possible rotations
    eps *= 0.5;
    eps += eps.transpose();

    // rotate to calculation system
    eps = _rotation * eps * _rotation.transpose();

    strain(0,0) = eps(1,1);
    strain(1,1) = eps(2,2);
    strain(2,2) = eps(3,3);
    strain(1,2) = strain(2,1) = eps(3,2);
    strain(0,2) = strain(2,0) = eps(3,1);
    strain(0,1) = strain(1,0) = eps(2,1);
    
  }

  return compatible;
}




void
BulkCrystal::print_info(void) const
{
  std::ostringstream os;
  os << "Bravais lattice : "
     << CrystalDefs::bravais_short_to_long_name(_bravais_lattice)
     << " (" << get_lattice_type() << ")" << std::endl;
  os << "a = " << _lattice_constant[0];
  if (_lattice_constant[1] != _lattice_constant[0])
    os << ", b = " << _lattice_constant[1];
  if (_lattice_constant[2] != _lattice_constant[1])
    os << ", c = " << _lattice_constant[2];
  os << std::endl;

  if ((_angles[0] != 90.0) || (_angles[1] != 90.0) || (_angles[2] != 90.0))
  {
    os << "alpha = " << _angles[0] << ", ";
    os << "beta = " << _angles[1] << ", ";
    os << "gamma = " << _angles[2] << "\n";
  }

  std::string mil = _options.get_option("x-growth-direction", " ");
  boost::replace_all(mil, ",", " ");
  os << "x growth direction : " << mil << std::endl;
  mil = _options.get_option("y-growth-direction", " ");
  boost::replace_all(mil, ",", " ");
  os << "y growth direction : " << mil << std::endl;
  mil = _options.get_option("z-growth-direction", " ");
  boost::replace_all(mil, ",", " ");
  os << "z growth direction : " << mil << std::endl;
  os << "Euler angles :" << " alpha = " << _euler_angles[0] / M_PI * 180 << ","
                         << " beta = " << _euler_angles[1] / M_PI * 180 << ","
                         << " gamma = " << _euler_angles[2] / M_PI * 180;
  Messages::info(os.str());
}

void
BulkCrystal::get_miller_indices(const Tensor1& vec, std::vector<int>& miller) const
{
  miller = std::vector<int>(3, 0);

  double h = vec * _conv_vec(1);
  double k = vec * _conv_vec(2);
  double l = vec * _conv_vec(3);

  double hh = abs(h);
  hh = (hh < 1e-3) ? 1.0 : hh;
  double kk = abs(k);
  kk = (kk < 1e-3) ? 1.0 : kk;
  double ll = abs(l);
  ll = (ll < 1e-3) ? 1.0 : ll;
  double minval = std::min(hh, std::min(kk, ll)); 

  h /= minval;
  k /= minval;
  l /= minval;

  double mul = 1000;
  
  double m0 = std::round(h * mul) / 10.0;
  double m1 = std::round(k * mul) / 10.0;
  double m2 = std::round(l * mul) / 10.0;

  miller[0] = std::trunc(m0);
  miller[1] = std::trunc(m1);
  miller[2] = std::trunc(m2);

  int div = std::gcd(miller[0], std::gcd(miller[1], miller[2]));
  miller[0] /= div;
  miller[1] /= div;
  miller[2] /= div;
  
}


void
BulkCrystal::calculate_euler_angles(void)
{
  // calculate Euler angles

  // formulas taken from en.wikipedia.org/wiki/Euler_angles
  // and https://eecs.qmul.ac.uk/~gslabaugh/publications/euler.pdf

  // Note: they are calculated from the transpose of _rotation
  double alpha = 0;
  double beta = acos(_rotation(3, 3));
  double gamma = 0;

  if (abs(_rotation(3,3)) < (1.0 - 1e-6))
  {

    double sb1 = sin(beta);
    alpha = atan2(-_rotation(3,2)/sb1, -_rotation(3,1)/sb1);


    if (abs(alpha) > M_PI_2)
    {
      double sb2 = -sb1;
      alpha = atan2(-_rotation(3,2)/sb2, -_rotation(3,1)/sb2);
      gamma = atan2(-_rotation(2,3)/sb2, _rotation(1,3)/sb2);
      beta = -beta;
    }
    else
      gamma = atan2(-_rotation(2,3)/sb1, _rotation(1,3)/sb1);
  }
  else // R33 = +/-1
  {
    // set gamma = 0 arbitrarily

    if (_rotation(3,3) < 0) // R33 = -1
    {
      alpha = atan2(-_rotation(2,1), _rotation(2,2));
    }
    else
    {
      alpha = atan2(_rotation(1,2), _rotation(1,1));
    }
  }

  _euler_angles[0] = alpha;
  _euler_angles[1] = beta;
  _euler_angles[2] = gamma;
}


void
BulkCrystal::get_euler_angles(double& alpha, double& beta, double& gamma) const
{
  alpha = _euler_angles[0];
  beta = _euler_angles[1];
  gamma = _euler_angles[2];
}
void



BulkCrystal::read_database(void)
{

  if ( !(_mat->is_alloy()) )
  {
    //lattice constant are expressed in Angstrom

    Database db = _mat->get_database();
    
    _lattice_type = db.get("structure", "none");

    db.set_section("lattice");
    _lattice_constant[0] = db.get("a", 0.0) * 10.0;
    //if (_lattice_constant[0] == 0.0) Messages::error("At least "
    //    "lattice constant a must be defined !!!!");
    _lattice_constant[1] = db.get("b", 0.0) * 10.0;
    if (_lattice_constant[1] == 0.0) _lattice_constant[1] = _lattice_constant[0];
    _lattice_constant[2] = db.get("c", 0.0) * 10.0;
    if (_lattice_constant[2] == 0.0) _lattice_constant[2] = _lattice_constant[0];
    _angles[0] = db.get("alpha", 90.0);
    _angles[1] = db.get("beta", 90.0);
    _angles[2] = db.get("gamma", 90.0);

    db.set_section("atomistic_structure");

    unsigned int n_basis_specie = db.get("n_basis_specie", 0);

    // the unique label for atoms in the primitive cell
    Atom::label_t label = 0;

    for (unsigned int i = 1; i <= n_basis_specie; i++)
    {
      std::string record, s, n_s;
      std::stringstream out;

      out << i;
      s = out.str();

      record = "n_" + s;
      unsigned int n = db.get(record.c_str(), 0);
      record = "specie_" + s;
      std::string specie = db.get(record.c_str(), "H");

      for (unsigned int j = 1; j <= n; j++)
      {
        record.clear(); n_s.clear();
        record = "T_" + s + "_";
        out.str(std::string());
        out.clear(std::stringstream::goodbit);
        out << j;
        n_s = out.str();
        n_s = record + n_s;

        //Putting specie label (defined by an integer) in label data
        //It's used in cut_and_change_specie() and build_random_alloy()
        label++;
        Atom tmp;
        tmp.set_label(label);
        tmp.set_specie(specie);

        libMesh::RealVectorValue v(0.0);
     
        // offers two alternative ways of setting vectors
        // either with component _a _b _c or as vectors 
        if (db.has_variable(n_s+"_a"))
        {
          v(0) = db.get(n_s+"_a", 0.0);
          v(1) = db.get(n_s+"_b", 0.0);
          v(2) = db.get(n_s+"_c", 0.0);
        }
        else
        {
          db.get(n_s, v);
        } 

        tmp.set_position(v);
        _lattice_basis.push_back(tmp);
      }
    }
  }

  if (_mat->is_alloy())
  {
    //Cannot act dynamic cast on mat itself because constant
    const Alloy* mat_alloy = dynamic_cast<const Alloy*>(_mat);
    double molar_fraction = _mat->get_options().get_option("x", 1.0);

    //NOTE: IMPLEMENTATION IS GOOD ONLY FOR BINARY COMPOUNDS
    //Nota: usiamo solo un pointer perche' per tutti i materiali viene istanziato solo un oggetto
    //database, a cui di volta in volta (ogni volta che chiamiamo un get) viene associato un data file.
    //Quindi non possiamo inizializzare 3 oggetti database e portarceli appresso, perche' saranno tutti
    //collegati al datafile settato dall'ultima assegnazione. Questa cosa va cambiata nella classe Database
    // (TODO)
    //
    // 2014-01-19 Matthias: I dont' think this is the case anymore??

    //This take two lines as I don't know how to specify that i call non constant method
    //if I do Database* db = &(mat_alloy->get_database())
    Database tmp_db = mat_alloy->get_database();
    Database* db = &tmp_db;

    if (db->get("alloy", 2) == 2)
    {
      double ax_1, ay_1, az_1, ax_2, ay_2, az_2;
      db = &(mat_alloy->get_component_A()->get_database());
      db->set_section("lattice");

      //We express lattice constant in Angstrom
      ax_1 = db->get("a", 0.0) * 10.0;
      if (ax_1 == 0.0) std::cerr << "At least lattice constant a must be defined !!!!" << std::endl;

      ay_1 = db->get("b", 0.0) * 10.0;
      if (ay_1 == 0.0) ay_1 = ax_1;

      az_1 = db->get("c", 0.0) * 10.0;
      if (az_1 == 0.0) az_1 = ax_1;

      db = &(mat_alloy->get_component_B()->get_database());
      db->set_section("lattice");

      ax_2 = db->get("a", 0.0) * 10.0;
      if (ax_2 == 0.0) std::cerr << "At least lattice constant a must be defined !!!!" << std::endl;

      ay_2 = db->get("b", 0.0) * 10.0;
      if (ay_2 == 0.0) ay_2 = ax_2;

      az_2 = db->get("c", 0.0) * 10.0;
      if (az_2 == 0.0) az_2 = ax_2;

      //Setting lattice parameters for the alloy
      _lattice_constant[0] = ax_1 * molar_fraction + ax_2 * (1.0 - molar_fraction);
      _lattice_constant[1] = ay_1 * molar_fraction + ay_2 * (1.0 - molar_fraction);
      _lattice_constant[2] = az_1 * molar_fraction + az_2 * (1.0 - molar_fraction);

    }
    else
      throw RuntimeException("Could not initialize bulk for non-binary alloy");
      

    db = &tmp_db;
    db->set_section("");

    db->set_section("atomistic_structure");

    unsigned int n_basis_specie = db->get("n_basis_specie", 0);

    // Read in basis vectors
    Database* dbB = &(mat_alloy->get_component_B()->get_database());
    Database* dbA = &(mat_alloy->get_component_A()->get_database());

    dbA->set_section("");
    _lattice_type = dbA->get("structure", "none");

    dbB->set_section("atomistic_structure");
    dbA->set_section("atomistic_structure");

    // the unique label for atoms in the primitive cell
    Atom::label_t label = 0;
        
    for (unsigned int i = 1; i <= n_basis_specie; i++)
    {
      std::string record("");
      std::string s("");
      std::stringstream out;
      out << i;
      s = out.str();
      record = "n_" + s;
      unsigned int n_x = (dbB->get(record, 0));
      if (dbA->get(record,0) != n_x) 
        Messages::error("Alloy bulks do not correspond");

      record = "specie_" + s;
      std::string specie = dbA->get(record.c_str(), "H");
      
      for (unsigned int j = 1; j <= n_x; j++)
      {
        std::string s2;
        record = "T_" + s + "_";
        out.str(std::string());
        out.clear(std::stringstream::goodbit);
        out << j;
        s2 = out.str();
        s2 = record + s2;

        //Putting specie label (defined by an integer) in label data
        //It's used in cut_and_change_specie() and build_random_alloy()
        label++;
        Atom tmp;
        tmp.set_label(label);
        tmp.set_specie(specie);

        libMesh::RealVectorValue Ta(0.0);
        libMesh::RealVectorValue Tb(0.0);
        record = s2 + "_a";
        Tb(0) = dbB->get(record, 0.0);
        Ta(0) = dbA->get(record, 0.0);
        record = s2 + "_b";
        Tb(1) = dbB->get(record, 0.0);
        Ta(1) = dbA->get(record, 0.0);
        record = s2 + "_c";
        Tb(2) = dbB->get(record, 0.0);
        Ta(2) = dbA->get(record, 0.0);

        Ta *= molar_fraction;
        Ta += Tb * (1.0 - molar_fraction);
        //libMesh::RealVectorValue T = Ta * (molar_fraction) + Tb * (1.0 - molar_fraction);
        tmp.set_position(Ta);
        _lattice_basis.push_back(tmp);

      }
    }

  } 

  _bravais_lattice = CrystalDefs::get_bravais_lattice(_lattice_type);

}


void
BulkCrystal::set_cell_vectors(void)
{

  Tensor2Gen prim_vec_dir(0);
  _prim_vec = 0.0;
  _conv_vec = 0.0;

  if (_bravais_lattice.at(0) == 'a')
  {
  }

  else if (_bravais_lattice.at(0) == 'm')
  {
    // following the definitions in Ponce et al. Phys. Rev. Research 2, 033102 (2020)
    _conv_vec(1,1) = _lattice_constant[0];
    _conv_vec(2,1) = 0.0;
    _conv_vec(3,1) = 0.0;

    _conv_vec(1,2) = 0.0;
    _conv_vec(2,2) = _lattice_constant[1];
    _conv_vec(3,2) = 0.0;

    double beta = M_PI * _angles[1] / 180;
    _conv_vec(1,3) = _lattice_constant[2] * cos(beta);
    _conv_vec(2,3) = 0.0;
    _conv_vec(3,3) = _lattice_constant[2] * sin(beta);
  }

  else if (_bravais_lattice.at(0) == 'o')
  {
    _conv_vec(1,1) = _lattice_constant[0];
    _conv_vec(2,2) = _lattice_constant[1];
    _conv_vec(3,3) = _lattice_constant[2];

    if (_bravais_lattice.at(1) == 'I')
    {
      Tensor1 prim_vec1 = 0.5 * (-_conv_vec(1) + _conv_vec(2) + _conv_vec(3));
      Tensor1 prim_vec2 = 0.5 * ( _conv_vec(1) - _conv_vec(2) + _conv_vec(3));
      Tensor1 prim_vec3 = 0.5 * ( _conv_vec(1) + _conv_vec(2) - _conv_vec(3));

      _prim_vec(1, 1) = prim_vec1(1); _prim_vec(2, 1) = prim_vec1(2); _prim_vec(3, 1) = prim_vec1(3);
      _prim_vec(1, 2) = prim_vec2(1); _prim_vec(2, 2) = prim_vec2(2); _prim_vec(3, 2) = prim_vec2(3);
      _prim_vec(1, 3) = prim_vec3(1); _prim_vec(2, 3) = prim_vec3(2); _prim_vec(3, 3) = prim_vec3(3);
    }
    else if (_bravais_lattice.at(1) == 'F')
    {
      Tensor1 prim_vec1 = 0.5 * (_conv_vec(2) + _conv_vec(3));
      Tensor1 prim_vec2 = 0.5 * (_conv_vec(1) + _conv_vec(3));
      Tensor1 prim_vec3 = 0.5 * (_conv_vec(1) + _conv_vec(2));

      _prim_vec(1, 1) = prim_vec1(1); _prim_vec(2, 1) = prim_vec1(2); _prim_vec(3, 1) = prim_vec1(3);
      _prim_vec(1, 2) = prim_vec2(1); _prim_vec(2, 2) = prim_vec2(2); _prim_vec(3, 2) = prim_vec2(3);
      _prim_vec(1, 3) = prim_vec3(1); _prim_vec(2, 3) = prim_vec3(2); _prim_vec(3, 3) = prim_vec3(3);
    }
  }

  else if (_bravais_lattice.at(0) == 'c')
  {
    _conv_vec = _lattice_constant[0] * Tensor2Gen(1.0);

    if (_bravais_lattice.at(1) == 'I')
    {
      prim_vec_dir(1, 1) = -0.5; prim_vec_dir(2, 1) = 0.5; prim_vec_dir(3, 1) = 0.5;
      prim_vec_dir(1, 2) = 0.5; prim_vec_dir(2, 2) = -0.5; prim_vec_dir(3, 2) = 0.5;
      prim_vec_dir(1, 3) = 0.5; prim_vec_dir(2, 3) = 0.5; prim_vec_dir(3, 3) = -0.5;

      _prim_vec = prim_vec_dir * _lattice_constant[0];
    }
    else if (_bravais_lattice.at(1) == 'F')
    {
      prim_vec_dir(1,1) = 0.0; prim_vec_dir(2,1) = 0.5; prim_vec_dir(3,1) = 0.5;
      prim_vec_dir(1,2) = 0.5; prim_vec_dir(2,2) = 0.0; prim_vec_dir(3,2) = 0.5;
      prim_vec_dir(1,3) = 0.5; prim_vec_dir(2,3) = 0.5; prim_vec_dir(3,3) = 0.0;

      _prim_vec = prim_vec_dir * _lattice_constant[0];
    }
  }

  else if (_bravais_lattice.compare("hP") == 0)
  {

    assert(_lattice_constant[0] == _lattice_constant[1]);

    prim_vec_dir(1,1) = 0.5; prim_vec_dir(2,1) = -sqrt(3.0) / 2.0; prim_vec_dir(3,1) = 0.0;
    prim_vec_dir(1,2) = 0.5; prim_vec_dir(2,2) = sqrt(3.0) / 2.0; prim_vec_dir(3,2) = 0.0;
    prim_vec_dir(1,3) = 0.0; prim_vec_dir(2,3) = 0.0; prim_vec_dir(3,3) = 1.0;

    _conv_vec(1,1) = prim_vec_dir(1,1) * _lattice_constant[0]; 
    _conv_vec(2,1) = prim_vec_dir(2,1) * _lattice_constant[0];
    _conv_vec(3,1) = prim_vec_dir(3,1) * _lattice_constant[0];
    _conv_vec(1,2) = prim_vec_dir(1,2) * _lattice_constant[0]; 
    _conv_vec(2,2) = prim_vec_dir(2,2) * _lattice_constant[0];
    _conv_vec(3,2) = prim_vec_dir(3,2) * _lattice_constant[0];
    _conv_vec(1,3) = prim_vec_dir(1,3) * _lattice_constant[2];
    _conv_vec(2,3) = prim_vec_dir(2,3) * _lattice_constant[2];
    _conv_vec(3,3) = prim_vec_dir(3,3) * _lattice_constant[2];

  }
  
  else if (_bravais_lattice.compare("hR") == 0)
  { 
  }

  else if (_bravais_lattice.at(0) == 't')
  {
    assert(_lattice_constant[0] == _lattice_constant[1]);

    _conv_vec(1,1) = _lattice_constant[0];
    _conv_vec(2,2) = _lattice_constant[0];
    _conv_vec(3,3) = _lattice_constant[2];

    if (_bravais_lattice.at(1) == 'I')
    {
      _prim_vec(1,1) = _conv_vec(1,1);
      _prim_vec(2,1) = _conv_vec(2,1);
      _prim_vec(3,1) = _conv_vec(3,1);
      _prim_vec(1,2) = _conv_vec(1,2);
      _prim_vec(2,2) = _conv_vec(2,2);
      _prim_vec(3,2) = _conv_vec(3,2);
      Tensor1 prim_vec = 0.5 * (_conv_vec(1) + _conv_vec(2) + _conv_vec(3));
      _prim_vec(1,3) = prim_vec(1);
      _prim_vec(2,3) = prim_vec(2);
      _prim_vec(3,3) = prim_vec(3);
    }
  }

  
  
  if (_bravais_lattice.at(1) == 'S')
  {
    Tensor1 prim_vec1 = 0.5 * (_conv_vec(1) - _conv_vec(2));
    Tensor1 prim_vec2 = 0.5 * (_conv_vec(1) + _conv_vec(2));
    Tensor1 prim_vec3 = _conv_vec(3);
    _prim_vec(1, 1) = prim_vec1(1); _prim_vec(2, 1) = prim_vec1(2); _prim_vec(3, 1) = prim_vec1(3);
    _prim_vec(1, 2) = prim_vec2(1); _prim_vec(2, 2) = prim_vec2(2); _prim_vec(3, 2) = prim_vec2(3);
    _prim_vec(1, 3) = prim_vec3(1); _prim_vec(2, 3) = prim_vec3(2); _prim_vec(3, 3) = prim_vec3(3);
  }
  else if (_bravais_lattice.at(1) == 'F')
  {

  }
  else if ((_bravais_lattice.at(1) == 'P') ||
           (_bravais_lattice.at(1) == 'R'))
  {
    _prim_vec = _conv_vec;
  }

  // create reciprocal lattice vectors
  Tensor1 a = _conv_vec(1);
  Tensor1 b = _conv_vec(2);
  Tensor1 c = _conv_vec(3);

  double vol = a * (b ^ c);

  _reciprocal_lattice[0] = (b ^ c)  / vol;
  _reciprocal_lattice[1] = (c ^ a)  / vol;
  _reciprocal_lattice[2] = (a ^ b)  / vol;

}

