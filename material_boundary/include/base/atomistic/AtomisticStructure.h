// $Id$

#ifndef _ATOMISTICSTRUCTURE_H_
#define _ATOMISTICSTRUCTURE_H_

#include "tensor.h"
#include "ModelOptions.h"
#include "TypeDefs.h"
#include "Control.h"
#include "Device.h"
#include "Database.h"
#include "Atom.h"
#include "BondMap.h"

//! Contains all needed data for an atomic structure
/*!
 *For any Atomistic region an atomistic structure is
 *built, containing atom positions and other informations
 *needed by any atomistic physical model.
 */
class AtomisticStructure
{
public:

//! A class for Dftb options
  class AtomisticStructureOptions
  {
  public:

    //!Options inizialization
    AtomisticStructureOptions(void);

    //! Copy constructor
    //AtomisticStructure(const AtomisticStructure &start);

    //!Options destructor
    ~AtomisticStructureOptions(void);

        //! Tells if structure has been passivated
    bool is_passivated;

    //!Tells if structure contains bond map informations
    bool contains_bond_map;

    //!Tells if structure has to be considered periodical;
    bool is_periodical;

    //!Tells if association has been already done
    bool is_associated;
  };


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


  //! Get set of regions covered by atomistic structure
  const std::set<std::string>& get_region(void);

  //! Get set of regions covered by atomistic structure (IDs)
  const std::set<ID>& get_IDset(void);

   //! Set the device we're working with
  void set_device(Device* const device);

  //! Get Device reference
  Device*  get_device(void);

  //! Return a const reference to structure atoms
  const std::vector<Atom>& get_structure_atoms(void) const;

  //! Return a writable reference to structure atoms
  std::vector<Atom>& get_structure_atoms(void);

  //! Get scale factor (from mesh_units to amstrong mesh_units/1e-10)
  const double& get_scale(void);

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

  //! Print upg file (etb_dataset may be changed into type options)
  void print_upg(const std::string& path, const std::string& etb_dataset);

  // get the number of atoms in the structure
  int get_N_atoms() const {return N_atoms;}

  // set number of atoms
  void set_N_atoms(long N) {N_atoms=N;}

  // get the number of types in the structure
  int get_N_types() const {return N_types;}

  // set number of atoms
  void set_N_types(int N) {N_types=N;}

  // say if it is periodic
  bool is_periodic() const {return _atomistic_structure_options.is_periodical;}

  //! Get atom type index (types are stored in _atom_types)
  int get_type_index   (  const std::string &   type    );

  void set_periodicity_vectors(const Tensor2Gen& T);

  void clear_atom_types(){_atom_types.clear();};

  void set_atom_types(const std::set<std::string>& types);

  //! Set Bond Map pointer
  void set_bondmap(BondMap* bondmap);

  //! Get bond map
  std::vector<std::vector<unsigned int> > const get_bond_map(void);

  //! AtomisticStructureOptions object pointer
  AtomisticStructureOptions _atomistic_structure_options;

  // Get element->atoms map
  //std::map<const Elem*, std::vector<unsigned int> >& get_elem_to_atoms(void);
  
  //! Get number of non hydrogen atoms
  unsigned int get_N_without_H(void);

private:


  //!Build element to atoms association map
  void build_elem_to_atoms(void);

  //! Associate at any alement atoms contained
      std::map<const Elem*, std::vector<unsigned int> > _elem_to_atoms;

  //!Associate elements: any atom keep tracks of the elements he belongs to
  void associate_elements();

  //! Scale factor (from mesh_units to amstrong mesh_units/1e-10)
  double _scale;

  //! Bond Map object pointer
  BondMap* _bondmap;

  //! Set the model options
  void set_options(const ModelOptions& options);

  //! Read structure from file
  void read_structure(const std::string& path);

  //! Read tgn structure
  void read_tgn(const std::string& path);

  //! Number of atoms in structure
  unsigned int N_atoms;

  //! Number of species
  int N_types;

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

  //! Reference material
  Material* _reference_material;

  //! Database of reference material
  Database _reference_material_db;

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
const double& AtomisticStructure::get_scale(void)
{
  return _scale;
}


inline
const std::vector<Atom>& AtomisticStructure::get_structure_atoms(void) const
{
  return _structure_atoms;
}


inline
std::vector<Atom>& AtomisticStructure::get_structure_atoms(void)
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
void
AtomisticStructure::set_bondmap(BondMap* bondmap)
{
  _bondmap = bondmap;
}

inline
std::vector<std::vector<unsigned int> > const
AtomisticStructure::get_bond_map()
{
  return _bondmap->get_bond_map();
}

//inline
//std::map<const Elem*, std::vector<unsigned int> >&
//AtomisticStructure::get_elem_to_atoms(void)
//{
//return _elem_to_atoms;
//}

#endif // _ATOMISTICSTRUCTURE_H_
