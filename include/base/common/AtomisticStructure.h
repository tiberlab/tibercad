#ifndef _ATOMISTICSTRUCTURE_H_
#define _ATOMISTICSTRUCTURE_H_

#include "tensor.h"
#include "ModelOptions.h"
#include "TypeDefs.h"
#include "Control.h"
#include "Device.h"
//#include "AtomisticGenerator.h"
//#include "AtomisticGenerator1D.h"
#include "Atom.h"

//C++ includes
#include <vector>
#include <set>
#include<iostream>
#include<fstream>
#include<sstream>
#include <map>



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

  ~AtomisticStructure(void);

  //! Get the structure options
  ModelOptions& get_options(void);

  //! Get the structure name
  const std::string& get_name(void);

  //! Get set of regions covered by atomistic structure
  const std::set<std::string>& get_region(void);

  //! Get set of regions covered by atomistic structure (IDs) 
  const std::set<ID>& get_IDset(void);

  //! Create a material with name /c name
  static AtomisticStructure* create(const std::string& name);

  //! Create a material with name /c name and options /c options
  static AtomisticStructure* create(const std::string& name, const ModelOptions& options);

  //! Set the device we're working with
  void set_device(Device* device);

  //! Get Device reference
  Device*  get_device(void);

  //! Return a reference to structure atoms
  const std::vector<Atom>& get_structure_atoms(void);

  //! Assign structure atoms vector
  void set_structure_atoms(const std::vector<Atom>& atoms);

  //! Return a reference to atom types
  const std::vector<std::string>& get_atom_types (void);

  //! Get periodicity vectors for the structure
  double* get_periodicity_vectors(void);

  //! Initialize the structure (up to now reading a structure from file is needed)
  void init(double a1 = 0.0, double a2 = 0.0, double a3 = 0.0); 

  //! Number of atoms in structure
  int N_atoms; 

  //! Number of species
  int N_types;

  //! Tell if atomistic structure has to be considered a periodical 
  //! structure (true) or a cluster (false)
  bool is_periodical;

  //! Get index of atom type
  int get_type_index(const std::string&);
  
  void add_atom_type(const std::string&);
  
  //! Vector containing structure atoms
  std::vector<Atom> _structure_atoms;

  //! Periodicity vectors in canonical basis
  double _periodicity_vectors[3][3];

  //! Set of all atom types in structure
  //std::vector<std::string> _atom_types;
  std::vector<std::string> _atom_types;


protected:

  //! Set the model options
  void set_options(const ModelOptions& options);

  //! Read structure from file
  void read_structure(const std::string& path);

  //! Print structure to file (format depends on extension used)
  void print_structure(const std::string& path);


private:

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
void AtomisticStructure::set_options(const ModelOptions& options)
{
  _options = options;

}

inline
ModelOptions& AtomisticStructure::get_options(void)
{
  return _options;
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
  return &_periodicity_vectors[0][0];
}


inline
Device* AtomisticStructure::get_device(void){
  return _device;
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


#endif // _ATOMISTICSTRUCTURE_H_
