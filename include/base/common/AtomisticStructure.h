

#ifndef _ATOMISTICSTRUCTURE_H_
#define _ATOMISTICSTRUCTURE_H_

#include "tensor.h"
#include "ModelOptions.h"

//C++ includes
#include <vector>
#include<iostream>
#include<fstream>
#include<sstream>


//! Contains atom definition
/*!
 *Atom is defined by atomic specie and a three component
 *vector (from library tensor.h) for the position.
 */
class atom
{
public:
  std::string specie;
  Tensor1 position;
};

//-------------------------------------------------
// Inline Memeber Functions
//-------------------------------------------------


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

  //! Create a material with name /c name
  static AtomisticStructure* create(const std::string& name);

  //! Create a material with name /c name and options /c options
  static AtomisticStructure* create(const std::string& name, const ModelOptions& options);

  //! Return a reference to structure atoms
  const std::vector<atom>& get_structure_atoms(void);

  //! Assign structure atoms vector
  void set_structure_atoms(const std::vector<atom>& atoms);

  //! Return a reference to atom types
  const std::vector<std::string>& get_atom_types (void);

  //! Initialize the structure (up to now reading a structure from file is needed)
  void init(void); 

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

  //! Vector containing structure atoms
  std::vector<atom> _structure_atoms;

  //! List of all atom types in structure
  std::vector<std::string> _atom_types;

  //! Number of atoms in structure
  unsigned int _N_atoms; 

  //! Tell if atomistic structure has to be considered a periodical 
  //! structure (true) or a cluster (false)
  bool _is_periodical;

  //! Periodicity vectors in canonical basis
  double _periodicity_vectors[3][3];

  //! Tell if the object has been already initialized
  bool _is_initialized;

};

//----------------------------------------------------
// Inline member functions
//----------------------------------------------------

inline
AtomisticStructure::AtomisticStructure(const std::string& name)
  :_name(name)
{ 
  // Default initializations
  _N_atoms = 0;
  _is_periodical = false;
  for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < 3; j++)
	{
	  _periodicity_vectors[i][j] = 0.0;
	}
    }

}


inline
AtomisticStructure::~AtomisticStructure(void)
{
}


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
const std::vector<atom>& AtomisticStructure::get_structure_atoms(void)
{
  return _structure_atoms;
}


inline
void AtomisticStructure::set_structure_atoms(const std::vector<atom>& atoms)
{
  _structure_atoms = atoms;
}


inline
const std::vector<std::string>& AtomisticStructure::get_atom_types(void)
{
  return _atom_types;
}


#endif // _ATOMISTICSTRUCTURE_H_
