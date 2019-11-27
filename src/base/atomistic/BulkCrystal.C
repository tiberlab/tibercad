#include "BulkCrystal.h"
#include "Database.h"
#include "Alloy.h"
#include "RotatedCrystal.h"
#include "RuntimeException.h"
#include <fstream>

BulkCrystal*
BulkCrystal::create(const Material* mat, const ModelOptions& options)
{
  return new BulkCrystal(mat, options);
}

BulkCrystal::BulkCrystal(const Material* mat, const ModelOptions& options)
:_lattice_constant(3, 0.0),
_angles(3, 90.0),  
_rotation(1),
_prim_vec(0),
_rotated_prim_vec(0)
{
  _options = ModelOptions(mat->get_options()); 
  _options += options;
  _mat = mat;
}

void
BulkCrystal::init(void)
{
 
  Atom tmp_atom;

  //Careful with the order of the calls, it's important
  read_database();
  set_prim_vec();

  _basis = _lattice_basis;
   for (std::vector<Atom>::iterator it = _basis.begin(); 
     it != _basis.end(); ++it)
   {
   it->set_position(_prim_vec * it->get_ttype_position());
   }

  //! Get rotation and apply to basis and primitive vectors
  build_rotation();
  set_ttype_lattice_vectors(_rotated_prim_vec);
  _atoms = _rotated_basis;

  refresh();
   
  
}


void
BulkCrystal::build_rotation(void)
{
  std::ostringstream os;

  Atom tmp_atom;

  //Retrieve rotation (now from material, in future this will 
  //be calculated directly here
  _rotation = _mat->get_rotated_crystal().RotMatrix; 
  
  //os << "Bulk Material " << _mat->get_name() << 
  //  " is created with a rotation " << std::endl<< _rotation << std::endl;
  //Messages::info(os.str());
  


  _rotated_prim_vec = _rotation * _prim_vec;

  //! Keep a copy of rotated crystal basis
   _rotated_basis = _basis;
   for (std::vector<Atom>::iterator it = _rotated_basis.begin(); 
     it != _rotated_basis.end(); ++it)
   {
   it->set_position(_rotation * it->get_ttype_position());
   }

}


void
BulkCrystal::read_database(void)
{
  Atom tmp;
  Tensor1 T;
  unsigned int j, n;
  std::string specie;

  if ( !(_mat->is_alloy()) )
  {
    //lattice constant are expressed in Amstrong

    Database db = _mat->get_database();
    db.set_section("lattice");
    _lattice_constant[0] = db.get("a", 0.0) * 10.0;
    if (_lattice_constant[0] == 0.0) Messages::error("At least "
        "lattice constant a must be defined !!!!");
    _lattice_constant[1] = db.get("b", 0.0) * 10.0;
    if (_lattice_constant[1] == 0.0) _lattice_constant[1] = _lattice_constant[0];
    _lattice_constant[2] = db.get("c", 0.0) * 10.0;
    if (_lattice_constant[2] == 0.0) _lattice_constant[2] = _lattice_constant[0];
    _angles[0] = db.get("alpha", 90.0);
    _angles[1] = db.get("beta", 90.0);
    _angles[2] = db.get("gamma", 90.0);

    db.set_section("atomistic_structure");

    _lattice_type = db.get("lattice_type", "none");

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
      n = db.get(record.c_str(), 0);
      record = "specie_" + s;
      specie = db.get(record.c_str(), "H");

      for (j = 1; j <= n; j++)
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
        tmp.set_label(label);
        tmp.set_specie(specie);

        std::vector<double> v(3,0.0);
     
        // offers two alternative ways of setting vectors
        // either with component _a _b _c or as vectors 
        if (db.has_variable(n_s+"_a"))
        {
          v[0] = db.get(n_s+"_a", 0.0);
          v[1] = db.get(n_s+"_b", 0.0);
          v[2] = db.get(n_s+"_c", 0.0);
        }
        else
        {
          db.get(n_s, v);
        } 
        T(1) = v[0]; 
        T(2) = v[1];
        T(3) = v[2];  
        tmp.set_position(T);
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

      //We express lattice constant in Amstrong
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
      throw RuntimeException("Could not initialize bulk for non-binar alloy");
      

    tmp_db = mat_alloy->get_database();
    db = &tmp_db;
    db->set_section("");
    db->set_section("atomistic_structure");
    _lattice_type = db->get("lattice_type", "none");

    tmp_db = mat_alloy->get_database();
    db = &tmp_db;

    db->set_section("atomistic_structure");
    unsigned int n_basis_specie = db->get("n_basis_specie", 0);

    // Read in basis vectors
    Tensor1 Tb, Ta;
    Database* dbB = &(mat_alloy->get_component_B()->get_database());
    Database* dbA = &(mat_alloy->get_component_A()->get_database());
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
      specie = dbA->get(record.c_str(), "H");
      
      for (j = 1; j <= n_x; j++)
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
        tmp.set_label(label);
        tmp.set_specie(specie);

        record = s2 + "_a";
        Tb(1) = dbB->get(record, 0.0);
        Ta(1) = dbA->get(record, 0.0);
        record = s2 + "_b";
        Tb(2) = dbB->get(record, 0.0);
        Ta(2) = dbA->get(record, 0.0);
        record = s2 + "_c";
        Tb(3) = dbB->get(record, 0.0);
        Ta(3) = dbA->get(record, 0.0);

        T = Ta * (molar_fraction) + Tb * (1.0 - molar_fraction);
        tmp.set_position(T);
        _lattice_basis.push_back(tmp);

      }
    }

  } 

}


void
BulkCrystal::set_prim_vec(void)
{

  Tensor2Gen prim_vec_dir(0);

  if (_lattice_type.compare("orthorombic") == 0) {

    prim_vec_dir(1,1) = 1.0; prim_vec_dir(2,1) = 0; prim_vec_dir(3,1) = 0;
    prim_vec_dir(1,2) = 0; prim_vec_dir(2,2) = 1.0; prim_vec_dir(3,2) = 0;
    prim_vec_dir(1,3) = 0; prim_vec_dir(2,3) = 0; prim_vec_dir(3,3) = 1.0;

    _prim_vec(1,1) = prim_vec_dir(1,1) * _lattice_constant[0];
    _prim_vec(2,2) = prim_vec_dir(2,2) * _lattice_constant[1];
    _prim_vec(3,3) = prim_vec_dir(3,3) * _lattice_constant[2];

  }

  if (_lattice_type.compare("tetragonal") == 0) {
 
    assert((_lattice_constant[0] == _lattice_constant[1]));

    prim_vec_dir(1,1) = 1.0; prim_vec_dir(2,1) = 0; prim_vec_dir(3,1) = 0;
    prim_vec_dir(1,2) = 0; prim_vec_dir(2,2) = 1.0; prim_vec_dir(3,2) = 0;
    prim_vec_dir(1,3) = 0; prim_vec_dir(2,3) = 0; prim_vec_dir(3,3) = 1.0;

    _prim_vec(1,1) = prim_vec_dir(1,1) * _lattice_constant[0];
    _prim_vec(2,2) = prim_vec_dir(2,2) * _lattice_constant[1];
    _prim_vec(3,3) = prim_vec_dir(3,3) * _lattice_constant[2];

  }

  if (_lattice_type.compare("cubic") == 0) {

    assert((_lattice_constant[0] == _lattice_constant[1]) && (_lattice_constant[1] == _lattice_constant[2]));

    prim_vec_dir(1,1) = 1.0 ; prim_vec_dir(2,1) = 0; prim_vec_dir(3,1) = 0;
    prim_vec_dir(1,2) = 0; prim_vec_dir(2,2) = 1.0; prim_vec_dir(3,2) = 0;
    prim_vec_dir(1,3) = 0; prim_vec_dir(2,3) = 0; prim_vec_dir(3,3) = 1.0;

    _prim_vec = prim_vec_dir * _lattice_constant[0];

  }

  else if (_lattice_type.compare("bcc") == 0) {

    assert((_lattice_constant[0] == _lattice_constant[1]) && (_lattice_constant[1] == _lattice_constant[2]));

    prim_vec_dir(1,1) = -0.5; prim_vec_dir(2,1) = 0.5; prim_vec_dir(3,1) = 0.5;
    prim_vec_dir(1,2) = 0.5; prim_vec_dir(2,2) = -0.5; prim_vec_dir(3,2) = 0.5;
    prim_vec_dir(1,3) = 0.5; prim_vec_dir(2,3) = 0.5; prim_vec_dir(3,3) = -0.5;

    _prim_vec = prim_vec_dir * _lattice_constant[0];

  }

  else if (_lattice_type.compare("fcc") == 0) {

    assert((_lattice_constant[0] == _lattice_constant[1]) && (_lattice_constant[1] == _lattice_constant[2]));

    prim_vec_dir(1,1) = 0.0; prim_vec_dir(2,1) = 0.5; prim_vec_dir(3,1) = 0.5;
    prim_vec_dir(1,2) = 0.5; prim_vec_dir(2,2) = 0.0; prim_vec_dir(3,2) = 0.5;
    prim_vec_dir(1,3) = 0.5; prim_vec_dir(2,3) = 0.5; prim_vec_dir(3,3) = 0.0;

    _prim_vec = prim_vec_dir * _lattice_constant[0];

  }

  else if ((_lattice_type.compare("hexagonal") == 0) ||
           (_lattice_type.compare("TMD") == 0)) {

    assert(_lattice_constant[0] == _lattice_constant[1]);

    prim_vec_dir(1,1) = 0.5; prim_vec_dir(2,1) = -sqrt(3.0) / 2.0; prim_vec_dir(3,1) = 0.0;
    prim_vec_dir(1,2) = 0.5; prim_vec_dir(2,2) = sqrt(3.0) / 2.0; prim_vec_dir(3,2) = 0.0;
    prim_vec_dir(1,3) = 0.0; prim_vec_dir(2,3) = 0.0; prim_vec_dir(3,3) = 1.0;

    _prim_vec(1,1) = prim_vec_dir(1,1) * _lattice_constant[0]; 
    _prim_vec(2,1) = prim_vec_dir(2,1) * _lattice_constant[0];
    _prim_vec(3,1) = prim_vec_dir(3,1) * _lattice_constant[0];
    _prim_vec(1,2) = prim_vec_dir(1,2) * _lattice_constant[0]; 
    _prim_vec(2,2) = prim_vec_dir(2,2) * _lattice_constant[0];
    _prim_vec(3,2) = prim_vec_dir(3,2) * _lattice_constant[0];
    _prim_vec(1,3) = prim_vec_dir(1,3) * _lattice_constant[2];
    _prim_vec(2,3) = prim_vec_dir(2,3) * _lattice_constant[2];
    _prim_vec(3,3) = prim_vec_dir(3,3) * _lattice_constant[2];

  }
  
  else if (_lattice_type.compare("anatase") == 0) {

    assert(_lattice_constant[0] == _lattice_constant[1]);

    prim_vec_dir(1,1) = 1.0; prim_vec_dir(2,1) = 0.0; prim_vec_dir(3,1) = 0.0;
    prim_vec_dir(1,2) = 0.0; prim_vec_dir(2,2) = 1.0; prim_vec_dir(3,2) = 0.0;
    prim_vec_dir(1,3) = 0.5; prim_vec_dir(2,3) = 0.5; prim_vec_dir(3,3) = 0.5;

    _prim_vec(1,1) = prim_vec_dir(1,1) * _lattice_constant[0];
    _prim_vec(2,2) = prim_vec_dir(2,2) * _lattice_constant[0];
    _prim_vec(1,3) = prim_vec_dir(1,3) * _lattice_constant[0]; 
    _prim_vec(2,3) = prim_vec_dir(2,3) * _lattice_constant[0]; 
    _prim_vec(3,3) = prim_vec_dir(3,3) * _lattice_constant[2];

  }
  else
  {
    Messages::error("Lattice type "+ _lattice_type + " doesn't exist in material "+_mat->get_name());
  }

  //Messages::info( "Bulk Material "+ _mat->get_name()+
  //                               " is created with primitive vectors ");
}

