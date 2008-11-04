#ifndef _ATOMISTICSTRUCTURE_H_
#define _ATOMISTICSTRUCTURE_H_

#include "tensor.h"
#include "ModelOptions.h"
#include "TypeDefs.h"
#include "Control.h"
#include "Device.h"
#include "Database.h"
#include "getpot.h"
#include "Atom.h"

//forward declaration
class BondMap;

//! Contains all needed data for an atomic structure
/*!
 *For any Atomistic region an atomistic structure is
 *built, containing atom positions and other informations
 *needed by any atomistic physical model.
 */
class AtomisticStructure
{
public:


  //! Constructor for AtomisticStructure class object
  /*!
   * A name to keep reference of structure must
   * be specified
   *
   */
  AtomisticStructure(const std::string& name);

  //! Destructor for AtomisticStructure class object
  ~AtomisticStructure(void);

 //! Create a material with name /c name
  static AtomisticStructure* create(const std::string& name);

  //! Create a material with name /c name and options /c options
  static AtomisticStructure* create(const std::string& name, const ModelOptions& options);

  //! Get the structure options
  ModelOptions& get_options(void);

  //! Get the structure name
  const std::string& get_name(void);

 //!Utility for parsing data from database (needed when Materials::init() has not been launched yet)
  GetPot getdata(const Material* mat);

  //! Get set of regions covered by atomistic structure
  const std::set<std::string>& get_region(void);

  //! Get set of regions covered by atomistic structure (IDs)
  const std::set<ID>& get_IDset(void);

   //! Set the device we're working with
  void set_device(Device* const device);

  //! Get Device reference
  Device*  get_device(void);

  //! Return a reference to structure atoms
  const std::vector<Atom>& get_structure_atoms(void);

  //! Assign structure atoms vector
  void set_structure_atoms(const std::vector<Atom>& atoms);

  //! Return a reference to atom types
  const std::vector<std::string>& get_atom_types (void);

  //! Get periodicity vectors for the structure:
  double* get_periodicity_vectors(void);

  //! Initialize the structure
  void init(void);

  //! Add a type to the atom types list
  void add_atom_type(const std::string&);

 //! Print structure to file (format depends on extension used)
  void print_structure(const std::string& path);

  //! A tool for printing atomic charges on output
  void print_structure(const std::string& path, double const* const charges);

  // get the number of atoms in the structure
  int get_N_atoms() const {return N_atoms;}

  // set number of atoms
  void set_N_atoms(long N) {N_atoms=N;}

  // get the number of types in the structure
  int get_N_types() const {return N_types;}

  // set number of atoms
  void set_N_types(int N) {N_types=N;}

  // say if it is periodic
  bool is_periodic() const {return is_periodical;}

  //! Get atom type index (types are stored in _atom_types)
  int get_type_index   (  const std::string &   type    );

  void set_periodicity_vectors(const Tensor2Gen& T);

  void clear_atom_types(){_atom_types.clear();};

  void set_atom_types(const std::set<std::string>& types);

  //! Set Bond Map pointer
  void set_bondmap(BondMap* bondmap);


private:

  //! Bond Map object pointer
  BondMap* _bondmap;

  //! Set the model options
  void set_options(const ModelOptions& options);

  //! Read structure from file
  void read_structure(const std::string& path);

  //! Number of atoms in structure
  unsigned int N_atoms;

  //! Number of species
  int N_types;

  //! Tell if atomistic structure has to be considered a periodical
  //! structure (true) or a cluster (false)
  bool is_periodical;

  //! Vector containing structure atoms
  std::vector<Atom> _structure_atoms;

  //! Periodicity vectors in canonical basis
  double _periodicity_vectors[9];

  //! Set of all atom types in structure
  std::vector<std::string> _atom_types;

  //! Options for the structure (from Atomistic Region)
  ModelOptions _options;

  //! Name of the structure (will be the same of associated Atomistic region)
  std::string _name;

  //! Set of regions covered by atomistic structure (names: TiberCAD
  //! defined regions)
  std::set <std::string> _regionset;

  //! Set of mesh regions covered by atomistic structure (numbers: mesh
  //! regiones)
  std::set <ID> _IDset;

  //! Tell if the object has been already initialized
  bool _is_initialized;

  //! Contains reference to device we're working with
  Device* _device;

};

//----------------------------------------------------
// Inline member functions
//----------------------------------------------------


inline
ModelOptions& AtomisticStructure::get_options(void)
{
  return _options;
}


inline
void AtomisticStructure::set_options(const ModelOptions& options)
{
  _options = options;

}


inline
const std::string& AtomisticStructure::get_name(void)
{
  return _name;
}


inline
const std::vector<Atom>& AtomisticStructure::get_structure_atoms(void)
{
  return _structure_atoms;
}


inline
void AtomisticStructure::set_structure_atoms(const std::vector<Atom>& atoms)
{
  _structure_atoms = atoms;
}


inline
const std::vector<std::string>& AtomisticStructure::get_atom_types(void)
{
  return _atom_types;
}


inline
double* AtomisticStructure::get_periodicity_vectors(void){
  return &_periodicity_vectors[0];
}


inline
Device* AtomisticStructure::get_device(void)
{
  return _device;
}

inline
void AtomisticStructure::set_device(Device* const device)
{
 _device = &(*device);
}


inline
const std::set<std::string>&
AtomisticStructure::get_region(void)
{
  return _regionset;
}


inline
const std::set<ID>&
AtomisticStructure::get_IDset(void)
{
  return _IDset;
}

inline
void AtomisticStructure::set_bondmap(BondMap* bondmap)
{
  _bondmap = bondmap;
}

#endif // _ATOMISTICSTRUCTURE_H_
