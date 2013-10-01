// $Id$

#ifndef _ATOMISTICSTRUCTURE_H_
#define _ATOMISTICSTRUCTURE_H_

#include "tensor.h"
#include "ModelOptions.h"
#include "TypeDefs.h"
#include "Device.h"
#include "Database.h"
#include "InitFailedException.h"
#include "AtomisticBasis.h"

class Material;
class UnstructuredMesh;

//! Contains all needed data for an atomic structure
/*!
 *For any Atomistic region an atomistic structure is
 *built, containing atom positions and other informations
 *needed by any atomistic physical model.
 */
class AtomisticStructure: public AtomisticBasis
{
public:

//! A class for Dftb options
  class AtomisticStructureOptions
  {
  public:

    //!Options inizialization
    AtomisticStructureOptions(void);

    //!Options destructor
    ~AtomisticStructureOptions(void);

    //!Tells if association has been already done
    bool is_associated;
  };



  //! Destructor for AtomisticStructure class object
  ~AtomisticStructure(void);

 //! Create a material with name /c name
  static AtomisticStructure* create();

  //! Get the structure options
  ModelOptions& get_options(void);

  //! Get the structure name
  const std::string& get_name(void);

  //! Get set of regions covered by atomistic structure
  const std::set<std::string>& get_region(void);

  //! Get set of regions covered by atomistic structure (IDs)
  const std::set<ID>& get_IDset(void) const;

   //! Set the device we're working with
  void set_device(const Device* const device);

  //! Get Device reference
  const Device*  get_device(void) const;

  //! Create a grid aligned on the atoms
  void create_conformal_grid(UnstructuredMesh& mesh) const;


  //! Get scale factor (from mesh_units to amstrong mesh_units/1e-10)
  const double& get_scale(void);

  //! Initialize a structure (to be read from input file)
  void init(const std::string& name, const Device* const device, const ModelOptions& options);

  // ! Initialize a structure (to be read from file)
  void init(const std::string& filename);

  //! Print structure to file (format depends on extension used)
  void print_structure(const std::string& path);

  //! A tool for printing atomic charges on output
  void print_structure(const std::string& path, double const* const charges);

  //! Print structure in internal tgn format (contains bondmap and element map)
  void print_tgn(const std::string& path) const;

  //! Print upg file (etb_dataset may be changed into type options)
  void print_upg(const std::string& path, const std::string& etb_dataset, 
                                          bool band_offsets = false);


  //! AtomisticStructureOptions object pointer
  AtomisticStructureOptions _atomistic_structure_options;

  // Get element->atoms map
  //std::map<const Elem*, std::vector<unsigned int> >& get_elem_to_atoms(void);
  
  //! Get number of non hydrogen atoms
  const unsigned int get_N_without_H(void) const;

  //! Set the model options
    void set_options(const ModelOptions& options);

  //! Get specie of atom i
  const Specie& get_specie(unsigned int i) const;

  //! Get atom Material
  const Material* get_material(const Atom& atom, bool parent = false) const;

  //! Get bond Material
  /*
   * If a couple of atoms belonging to a binary compound is specified,
   * the material can be assigned without ambiguity to the Anion.
   * The check is now hardcoded for a restricted class of materials (Nitrides
   * and Arsenides). With additional information in Material class, this can be made
   * more generic.
   */
  const Material* get_material(const Atom& atom1, const Atom& atom2,
      bool parent = false) const;

  //! Tells if random alloy approximation is used
  bool is_random_alloy(void);

  //! Tells if clustering should be applied
  bool build_clusters(void);

  //! Get the reference material
  const Material* get_reference_material(void) const;

  //! Apply reordering to atoms
  void reorder(const std::vector<unsigned int>& P);

private:

  //! Constructor for AtomisticStructure class object
  /*!
   * A name to keep reference of structure must
   * be specified
   *
   */
  AtomisticStructure();

  AtomisticStructure(const std::string& name);
  
  //! copy constructor
  AtomisticStructure(const AtomisticStructure& other);
  
  //! Calculate the number of atoms excluding hydrogens, usefule for 
    //! passivated structures
    void compute_N_without_H(void);

    //!Override lattice vectors from structure generation
    void parse_lattice_vectors(void); 

     //! Initialize the structure using mesh infos
    void init_mesh_structure(void);

    //Build mesh regions infos
    void parse_regions(void);

  //!Build element to atoms association map
  void build_elem_to_atoms(void);

  //! Associate at any alement atoms contained
      std::map<const Elem*, std::vector<unsigned int> > _elem_to_atoms;

  //!Associate elements: any atom keep tracks of the elements he belongs to
  void associate_elements();

  //! Scale factor (from mesh_units to amstrong mesh_units/1e-10)
  double _scale;

  //! Read structure from file
  void read_structure(const std::string& path);

  //! Read tgn structure
  void read_tgn(const std::string& path);


  //! Options for the structure (from Atomistic Region)
  ModelOptions _options;

  //! Name of the structure (will be the same of associated Atomistic region)
  std::string _name;

  //! Set of regions covered by atomistic structure (names: TiberCAD
  //! defined regions)
  //std::set <std::string> _regionset;

  //! Set of mesh regions covered by atomistic structure (numbers: mesh
  //! regiones)
  std::set <ID> _IDset;

  //! Tell if the object has been already initialized
  bool _is_initialized;

  //! True if the alloy must be treated as a Random Alloy
  //! (Species are used to distinguish the parent material the atom
  //! belongs to). By default VCA is used
  bool _random_alloy;

  //! True if clustering should be applied
  bool _clustering;

  //! Contains reference to device we're working with
  const Device* _device;

  //! Reference material
  Material* _reference_material;

  //! Database of reference material
  Database _reference_material_db;

  //!Number of atoms excluding hydrogens
  unsigned int _N_without_H;

  //! Manage structure printing
  void print_driver(void);

};

//----------------------------------------------------
// Inline member functions
//----------------------------------------------------


inline
ModelOptions& AtomisticStructure::get_options(void)
{
    assert( !(_options.is_empty()) );

    return _options;

 }

inline
bool AtomisticStructure::is_random_alloy(void)
{
  return _random_alloy;
}


inline
bool AtomisticStructure::build_clusters(void)
{
  return _clustering;
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
const Device* AtomisticStructure::get_device(void) const
{

  assert ( (_device != NULL) );

  return _device;

}

inline
void AtomisticStructure::set_device(const Device* const device)
{
 _device = device;
}


/*
inline
const std::set<std::string>&
AtomisticStructure::get_region(void)
{
  return _regionset;
}
*/


inline
const std::set<ID>&
AtomisticStructure::get_IDset(void) const
{
  return _IDset;
}


inline
const Specie&
AtomisticStructure::get_specie(unsigned int i) const
{
  return _atoms[i].get_specie();
}

inline
const unsigned int
AtomisticStructure::get_N_without_H(void) const
{
  return _N_without_H;
}


inline 
const Material* 
AtomisticStructure::get_reference_material(void) const
{
  return _reference_material;
}

//inline
//std::map<const Elem*, std::vector<unsigned int> >&
//AtomisticStructure::get_elem_to_atoms(void)
//{
//return _elem_to_atoms;
//}

#endif // _ATOMISTICSTRUCTURE_H_
